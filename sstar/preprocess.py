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
import scipy
import yaml
from multiprocessing import Process, Queue
from sstar.stats import *
from sstar.utils import read_data, filter_data, create_windows, multiprocessing_manager


def process_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, feature_config, is_phased, 
                 ploidy, output_dir, output_prefix, win_len, win_step, thread):
    """
    Description:
        Processes genotype data.

    Arguments:
        vcf_file str: Name of the VCF file containing genotype data.
        ref_ind_file str: Name of the file containing sample information from the reference population.
        tgt_ind_file str: Name of the file containing sample information from the target population.
        anc_allele_file str: Name of the file containing ancestral allele information.
        feature_config str: Name of the YAML file specifying what features should be used. 
        is_phased bool: True, genomes are phased; False, genomes are unphased.
        ploidy int: Ploidy of genomes.
        output_dir str: Directory storing the output files.
        output_prefix str: Prefix of the output files.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        thread int: Number of threads.
    """
    if feature_config is None:
        features = None
        output_genotypes = True
    else:
        with open(feature_config, 'r') as f:
            features = yaml.safe_load(f)

        features = features['features']
        output_genotypes = False

    ref_data, ref_samples, tgt_data, tgt_samples = read_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, is_phased)

    for c in tgt_data.keys():
        windows = create_windows(tgt_data[c]['POS'], c, win_step, win_len)

    res = multiprocessing_manager(worker_func=preprocess_worker, nrep=len(windows), thread=thread, windows=windows, ref_data=ref_data, tgt_data=tgt_data, 
                                  features=features, is_phased=is_phased, ploidy=ploidy, output_genotypes=output_genotypes)

    # x[0]: the chromosome name in number
    # x[1]: the start of the window
    # x[2]: the end of the window
    res.sort(key=lambda x: (x[0], x[1], x[2]))

    header = _create_header(ref_samples, tgt_samples, features, is_phased, ploidy, output_genotypes)
    _output(res, tgt_samples, header, features, is_phased, ploidy, output_dir, output_prefix, output_genotypes)


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

        if kwargs['output_genotypes'] is True:
            items['ref_gts'] = ref_gts
            items['tgt_gts'] = tgt_gts
        else:
            if 'number of private mutations' in kwargs['features'].keys():
                pvt_mut_nums = cal_pvt_mut_num(sub_ref_gts, sub_tgt_gts)
                items['pvt_mut_nums'] = pvt_mut_nums
            if 'individual allele frequency spectra' in kwargs['features'].keys():
                if kwargs['is_phased'] is True: ploidy = 1
                else: ploidy = kwargs['ploidy']
                spectra = cal_n_ton(tgt_gts, ploidy=ploidy)
                items['spectra'] = spectra
            if 'reference distances' in kwargs['features'].keys():
                ref_dists = cal_dist(ref_gts, tgt_gts)
                if 'all' in kwargs['features']['reference distances'].keys(): items['ref_dists'] = ref_dists
                if 'minimum' in kwargs['features']['reference distances'].keys(): items['min_ref_dists'] = np.min(ref_dists, axis=1)
                if 'maximum' in kwargs['features']['reference distances'].keys(): items['max_ref_dists'] = np.max(ref_dists, axis=1)
                if 'mean' in kwargs['features']['reference distances'].keys(): items['mean_ref_dists'] = np.mean(ref_dists, axis=1)
                if 'median' in kwargs['features']['reference distances'].keys(): items['median_ref_dists'] = np.median(ref_dists, axis=1)
                if 'variance' in kwargs['features']['reference distances'].keys(): items['var_ref_dists'] = np.var(ref_dists, axis=1)
                if 'skew' in kwargs['features']['reference distances'].keys(): items['skew_ref_dists'] = scipy.stats.skew(ref_dists, axis=1)
                if 'kurtosis' in kwargs['features']['reference distances'].keys(): items['kurtosis_ref_dists'] = scipy.stats.kurtosis(ref_dists, axis=1)
            if 'target distances' in kwargs['features'].keys():
                tgt_dists = cal_dist(tgt_gts, tgt_gts)
                if 'all' in kwargs['features']['target distances'].keys(): items['tgt_dists'] = tgt_dists
                if 'minimum' in kwargs['features']['target distances'].keys(): items['min_tgt_dists'] = np.min(tgt_dists, axis=1)
                if 'maximum' in kwargs['features']['target distances'].keys(): items['max_tgt_dists'] = np.max(tgt_dists, axis=1)
                if 'mean' in kwargs['features']['target distances'].keys(): items['mean_tgt_dists'] = np.mean(tgt_dists, axis=1)
                if 'median' in kwargs['features']['target distances'].keys(): items['median_tgt_dists'] = np.median(tgt_dists, axis=1)
                if 'variance' in kwargs['features']['target distances'].keys(): items['var_tgt_dists'] = np.var(tgt_dists, axis=1)
                if 'skew' in kwargs['features']['target distances'].keys(): items['skew_tgt_dists'] = scipy.stats.skew(tgt_dists, axis=1)
                if 'kurtosis' in kwargs['features']['target distances'].keys(): items['kurtosis_tgt_dists'] = scipy.stats.kurtosis(tgt_dists, axis=1)
            if 'sstar' in kwargs['features'].keys():
                sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 
                                                                     method=kwargs['features']['sstar']['genotype distance'], 
                                                                     match_bonus=kwargs['features']['sstar']['match bonus'],
                                                                     max_mismatch=kwargs['features']['sstar']['max mismatch'],
                                                                     mismatch_penalty=kwargs['features']['sstar']['mismatch penalty'])
                items['sstar'] = sstar_scores

        out_queue.put((chr_name, start, end, items))


def _create_header(ref_samples, tgt_samples, features, is_phased, ploidy, output_genotypes):
    """
    """
    if output_genotypes is True:
        if is_phased is True:
            haps = []
            for s in ref_samples:
                for i in range(ploidy):
                    haps.append(f'{s}_{i+1}')
            for s in tgt_samples:
                for i in range(ploidy):
                    haps.append(f'{s}_{i+1}')
            header = "\t".join(haps)
        else: header = "\t".join(ref_samples) + "\t" + "\t".join(tgt_samples)
    else:
        header = "chrom\tstart\tend\tsample"
        if 'sstar' in features.keys(): header += "\tS*_score"
        if 'number of private mutations' in features.keys(): header += "\tprivate_SNP_num"
        if 'individual allele frequency spectra' in features.keys():
            nsample = len(tgt_samples)*ploidy
            for i in range(nsample+1): header += f'\t{i}_ton'
        if 'reference distances' in features.keys():
            if 'minimum' in features['reference distances'].keys(): header += "\tmin_ref_dist"
            if 'maximum' in features['reference distances'].keys(): header += "\tmax_ref_dist"
            if 'mean' in features['reference distances'].keys(): header += "\tmean_ref_dist"
            if 'median' in features['reference distances'].keys(): header += "\tmedian_ref_dist"
            if 'variance' in features['reference distances'].keys(): header += "\tvar_ref_dist"
            if 'skew' in features['reference distances'].keys(): header += "\tskew_ref_dist"
            if 'kurtosis' in features['reference distances'].keys(): header += "\tkurtosis_ref_dist"
            if 'all' in features['reference distances'].keys(): 
                if is_phased is True: 
                    for s in ref_samples:
                        for i in range(ploidy):
                            header += f'\tref_dist_{s}_{i+1}'
                else:
                    for s in ref_samples:
                        header += f'\tref_dist_{s}'
        if 'target distances' in features.keys():
            if 'minimum' in features['target distances'].keys(): header += "\tmin_tgt_dist"
            if 'maximum' in features['target distances'].keys(): header += "\tmax_tgt_dist"
            if 'mean' in features['target distances'].keys(): header += "\tmean_tgt_dist"
            if 'median' in features['target distances'].keys(): header += "\tmedian_tgt_dist"
            if 'variance' in features['target distances'].keys(): header += "\tvar_tgt_dist"
            if 'skew' in features['target distances'].keys(): header += "\tskew_tgt_dist"
            if 'kurtosis' in features['target distances'].keys(): header += "\tkurtosis_tgt_dist"
            if 'all' in features['target distances'].keys(): 
                if is_phased is True: 
                    for s in tgt_samples:
                        for i in range(ploidy):
                            header += f'\ttgt_dist_{s}_{i+1}'
                else: 
                    for s in tgt_samples:
                        header += f'\ttgt_dist_{s}'

    return header


def _output(res, tgt_samples, header, features, is_phased, ploidy, output_dir, output_prefix, output_genotypes):
    """
    """
    os.makedirs(output_dir, exist_ok=True)
    if output_genotypes is True:
        os.makedirs(f'{output_dir}/genotypes', exist_ok=True)
        for r in res:
            chrom = r[0]
            start = r[1]
            end = r[2]
            items = r[3]
            os.makedirs(f'{output_dir}/genotypes/{chrom}', exist_ok=True)
            if is_phased is True: output_file = f'{output_dir}/genotypes/{chrom}/{output_prefix}.{chrom}.{start}-{end}.phased.genotypes'
            else: output_file = f'{output_dir}/genotypes/{chrom}/{output_prefix}.{chrom}.{start}-{end}.unphased.genotypes'
            with open(output_file, 'w') as f:
                f.write(f'{header}\n')
                for i in range(len(items['ref_gts'])):
                    ref_gts = "\t".join(items['ref_gts'][i].astype(str))
                    tgt_gts = "\t".join(items['tgt_gts'][i].astype(str))
                    f.write(f'{ref_gts}\t{tgt_gts}\n')
    else:
        output_file = f'{output_dir}/{output_prefix}.features'
        if is_phased is not True: ploidy = 1
        with open(output_file, 'w') as f:
            f.write(f'{header}\n')
            for r in res:
                chrom = r[0]
                start = r[1]
                end = r[2]
                items = r[3]
                for i in range(len(tgt_samples)*ploidy):
                    if ploidy != 1: sample = f'{tgt_samples[int(i/ploidy)]}_{i%ploidy+1}'
                    else: sample = tgt_samples[i]
                    out = ''
                    if 'sstar' in features.keys(): out += f'{items["sstar"][i]}'
                    if 'number of private mutations' in features.keys(): out += f'\t{items["pvt_mut_nums"][i]}'
                    if 'individual allele frequency spectra' in features.keys():
                        spectra = "\t".join(items["spectra"][i].astype(str))
                        out += f'\t{spectra}'
                    if 'reference distances' in features.keys():
                        if 'minimum' in features['reference distances'].keys(): out += f'\t{items["min_ref_dists"][i]}'
                        if 'maximum' in features['reference distances'].keys(): out += f'\t{items["max_ref_dists"][i]}'
                        if 'mean' in features['reference distances'].keys(): out += f'\t{items["mean_ref_dists"][i]}'
                        if 'median' in features['reference distances'].keys(): out += f'\t{items["median_ref_dists"][i]}'
                        if 'variance' in features['reference distances'].keys(): out += f'\t{items["var_ref_dists"][i]}'
                        if 'skew' in features['reference distances'].keys(): out += f'\t{items["skew_ref_dists"][i]}'
                        if 'kurtosis' in features['reference distances'].keys(): out += f'\t{items["kurtosis_ref_dists"][i]}'
                        if 'all' in features['reference distances'].keys(): 
                            dists = "\t".join(items["ref_dists"][i].astype(str))
                            out += f'\t{dists}'
                    if 'target distances' in features.keys():
                        if 'minimum' in features['target distances'].keys(): out += f'\t{items["min_tgt_dists"][i]}'
                        if 'maximum' in features['target distances'].keys(): out += f'\t{items["max_tgt_dists"][i]}'
                        if 'mean' in features['target distances'].keys(): out += f'\t{items["mean_tgt_dists"][i]}'
                        if 'median' in features['target distances'].keys(): out += f'\t{items["median_tgt_dists"][i]}'
                        if 'variance' in features['target distances'].keys(): out += f'\t{items["var_tgt_dists"][i]}'
                        if 'skew' in features['target distances'].keys(): out += f'\t{items["skew_tgt_dists"][i]}'
                        if 'kurtosis' in features['target distances'].keys(): out += f'\t{items["kurtosis_tgt_dists"][i]}'
                        if 'all' in features['target distances'].keys(): 
                            dists = "\t".join(items["tgt_dists"][i].astype(str))
                            out += f'\t{dists}'
                    f.write(f'{chrom}\t{start}\t{end}\t{sample}\t{out}\n')


if __name__ == '__main__':
    process_data(vcf_file="examples/data/real_data/sstar.example.biallelic.snps.vcf.gz", ref_ind_file="examples/data/ind_list/ref.ind.list", tgt_ind_file="examples/data/ind_list/tgt.ind.list", anc_allele_file=None, feature_config="examples/pre-trained/test.features.yaml", is_phased=True, ploidy=2, output_dir="sstar/test", output_prefix="test", win_len=50000, win_step=10000, thread=1)
    #process_data(vcf_file="./sstar/test/0/test.0.vcf", ref_ind_file="./sstar/test/0/test.0.ref.ind.list", tgt_ind_file="./sstar/test/0/test.0.tgt.ind.list", anc_allele_file=None, feature_config="examples/pre-trained/test.features.yaml", is_phased=True, ploidy=2, output_genotypes=False, output_dir="sstar/test", output_prefix="test", win_len=50000, win_step=10000, thread=1)
