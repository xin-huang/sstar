# Apache License Version 2.0
# Copyright 2023 Xin Huang
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
import numpy as np
import yaml
from multiprocessing import Process, Queue
from sstar.stats import *
from sstar.utils import read_data, filter_data, create_windows, multiprocessing_manager


def process_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, feature_file, output_dir, output_prefix, win_len, win_step, thread):
    """
    Description:
        Processes genotype data.

    Arguments:
        vcf_file str: Name of the VCF file containing genotype data.
        ref_ind_file str: Name of the file containing sample information from the reference population.
        tgt_ind_file str: Name of the file containing sample information from the target population.
        anc_allele_file str: Name of the file containing ancestral allele information.
        feature_file str: Name of the YAML file specifying what features should be used. 
        output_dir str: Directory storing the output files.
        output_prefix str: Prefix of the output files.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        thread int: Number of threads.
    """
    with open(feature_file, 'r') as f:
        features = yaml.safe_load(f)

    features = features['features']

    ref_data, ref_samples, tgt_data, tgt_samples = read_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, features['genotypes']['phased'])

    for c in tgt_data.keys():
        windows = create_windows(tgt_data[c]['POS'], c, win_step, win_len)

    res = multiprocessing_manager(worker_func=preprocess_worker, nrep=len(windows), thread=thread, windows=windows, ref_data=ref_data, tgt_data=tgt_data, features=features)

    # x[0]: the chromosome name in number
    # x[1]: the start of the window
    # x[2]: the end of the window
    res.sort(key=lambda x: (x[0], x[1], x[2]))

    if features['genotypes']['output']:
        header = _create_header(ref_samples, tgt_samples, features, features['genotypes']['output'])
        _output(res, header, features, output_dir, output_prefix, features['genotypes']['output'])


def preprocess_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        chr_name, start, end, ref_gts, tgt_gts, pos = in_queue.get()
        variants_not_in_ref = np.sum(ref_gts, axis=1) == 0
        sub_ref_gts = ref_gts[variants_not_in_ref]
        sub_tgt_gts = tgt_gts[variants_not_in_ref]
        sub_pos = pos[variants_not_in_ref]

        items = dict()

        if ('genotypes' in kwargs['features'].keys()) and (kwargs['features']['genotypes']['output'] is True):
            items['ref_gts'] = ref_gts
            items['tgt_gts'] = tgt_gts
        if ('number of private mutations' in kwargs['features'].keys()) and (kwargs['features']['number of private mutations']['output'] is True):
            pvt_mut_nums = cal_pvt_mut_num(sub_ref_gts, sub_tgt_gts)
            items['pvt_mut_nums'] = pvt_mut_nums
        if ('individual allele frequency spectra' in kwargs['features'].keys()) and (kwargs['features']['individual allele frequency spectra']['output'] is True):
            spectra = cal_n_ton(tgt_gts)
            items['spectra'] = spectra
        if 'pairwise distances' in kwargs['features'].keys():
            if ('reference and target populations' in kwargs['features']['pairwise distances'].keys()) and (kwargs['features']['pairwise distances']['reference and target populations']['output'] is True):
                ref_dists = cal_dist(ref_gts, tgt_gts)
                if kwargs['features']['pairwise distances']['reference and target populations']['all'] is True: items['all_ref_dists'] = ref_dists
                if kwargs['features']['pairwise distances']['reference and target populations']['minimum'] is True: items['min_ref_dists'] = np.min(ref_dists, axis=1)
                if kwargs['features']['pairwise distances']['reference and target populations']['maximum'] is True: items['max_ref_dists'] = np.max(ref_dists, axis=1)
                if kwargs['features']['pairwise distances']['reference and target populations']['mean'] is True: items['mean_ref_dists'] = np.mean(ref_dists, axis=1)
                if kwargs['features']['pairwise distances']['reference and target populations']['median'] is True: items['median_ref_dists'] = np.median(ref_dists, axis=1)
                if kwargs['features']['pairwise distances']['reference and target populations']['variance'] is True: items['var_ref_dists'] = np.var(ref_dists, axis=1)
                if kwargs['features']['pairwise distances']['reference and target populations']['skew'] is True: items['skew_ref_dists'] = np.skew(ref_dists, axis=1)
                if kwargs['features']['pairwise distances']['reference and target populations']['kurtosis'] is True: items['kurtosis_ref_dists'] = np.kurtosis(ref_dists, axis=1)
            if ('target population only' in kwargs['features'].keys()) and (kwargs['features']['pairwise distances']['target population only']['output'] is True):
                tgt_dists = cal_dist(tgt_gts, tgt_gts)
                if kwargs['features']['pairwise distances']['target population only']['all'] is True: items['all_tgt_dists'] = tgt_dists
                if kwargs['features']['pairwise distances']['target population only']['minimum'] is True: items['min_tgt_dists'] = np.min(tgt_dists, axis=1)
                if kwargs['features']['pairwise distances']['target population only']['maximum'] is True: items['max_tgt_dists'] = np.max(tgt_dists, axis=1)
                if kwargs['features']['pairwise distances']['target population only']['mean'] is True: items['mean_tgt_dists'] = np.mean(tgt_dists, axis=1)
                if kwargs['features']['pairwise distances']['target population only']['median'] is True: items['median_tgt_dists'] = np.median(tgt_dists, axis=1)
                if kwargs['features']['pairwise distances']['target population only']['variance'] is True: items['var_tgt_dists'] = np.var(tgt_dists, axis=1)
                if kwargs['features']['pairwise distances']['target population only']['skew'] is True: items['skew_tgt_dists'] = np.skew(tgt_dists, axis=1)
                if kwargs['features']['pairwise distances']['target population only']['kurtosis'] is True: items['kurtosis_tgt_dists'] = np.kurtosis(tgt_dists, axis=1)
        if ('sstar' in kwargs['features'].keys()) and (kwargs['features']['sstar']['output'] is True):
            sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 
                                                                 method=kwargs['features']['sstar']['genotype distance'], 
                                                                 match_bonus=kwargs['features']['sstar']['match bonus'],
                                                                 max_mismatch=kwargs['features']['sstar']['max mismatch'],
                                                                 mismatch_penalty=kwargs['features']['sstar']['mismatch penalty'])
            items['sstar'] = sstar_scores

        out_queue.put((chr_name, start, end, items))


def _create_header(ref_samples, tgt_samples, features, output_genotypes):
    """
    """
    if output_genotypes:
        if features['genotypes']['phased'] is True:
            haps = []
            for s in ref_samples:
                haps.append(s+'_hap1')
                haps.append(s+'_hap2')
            for s in tgt_samples:
                haps.append(s+'_hap1')
                haps.append(s+'_hap2')
            header = "\t".join(haps)
        else: header = "\t".join(ref_samples) + "\t" + "\t".join(tgt_samples)
    else:
        header = "chrom\tstart\tend"

    return header


def _output(res, header, features, output_dir, output_prefix, output_genotypes):
    """
    """
    os.makedirs(output_dir, exist_ok=True)
    if output_genotypes:
        os.makedirs(f'{output_dir}/genotypes', exist_ok=True)
        for r in res:
            chrom = r[0]
            start = r[1]
            end = r[2]
            items = r[3]
            if features['genotypes']['phased'] is True: output_file = f'{output_dir}/genotypes/{output_prefix}.{chrom}.{start}-{end}.phased.genotypes'
            else: output_file = f'{output_dir}/genotypes/{output_prefix}.{chrom}.{start}-{end}.unphased.genotypes'
            with open(output_file, 'w') as f:
                f.write(f'{header}\n')
                for i in range(len(items['ref_gts'])):
                    ref_gts = "\t".join(items['ref_gts'][i].astype(str))
                    tgt_gts = "\t".join(items['tgt_gts'][i].astype(str))
                    f.write(f'{ref_gts}\t{tgt_gts}\n')


if __name__ == '__main__':
    process_data(vcf_file="examples/data/real_data/sstar.example.biallelic.snps.vcf.gz", ref_ind_file="examples/data/ind_list/ref.ind.list", tgt_ind_file="examples/data/ind_list/tgt.ind.list", anc_allele_file=None, feature_file="examples/pre-trained/test.features.yaml", output_dir="sstar/test", output_prefix="test", win_len=50000, win_step=10000, thread=1)
