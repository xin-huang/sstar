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


def process_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, feature_file, output, win_len, win_step, thread):
    """
    Description:
        Processes genotype data.

    Arguments:
        vcf_file str: Name of the VCF file containing genotype data.
        ref_ind_file str: Name of the file containing sample information from the reference population.
        tgt_ind_file str: Name of the file containing sample information from the target population.
        anc_allele_file str: Name of the file containing ancestral allele information.
        feature_file str: Name of the YAML file specifying what features should be used. 
        output str: Name of the output file.
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

    print(res)


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


def _process(win_step, win_len, output, thread, **kwargs):
    """
    Description:
        Processes genotype data with the S* statistic and others.

    Arguments:
        win_step int: Length of sliding windows.
        win_len int: Step size of sliding windows.
        output str: Name of the output file.
        thread int: Number of threads.

    Keyword arguments:
        ref_data dict: Dictionary containing data from the reference population.
        tgt_data dict: Dictionary containing data from the target population.
        samples list: List containing sample information from the target population.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
    """
    ind_num = len(kwargs['samples'])
    header = 'chrom\tstart\tend\tsample\thap\t'
    header += '\t'.join([str(x)+'-ton' for x in range(ind_num*2+1)]) + '\t'
    header += '\t'.join(['pairwised_dist'+str(x+1) for x in range(ind_num*2)])
    header += '\tmean_pairwised_dist\tvar_pairwised_dist\tskew_pairwised_dist\tkurtosis_pairwised_dist'
    header += '\tmin_dist_to_ref\tS*_score\tprivate_SNP_num'
    worker_func = _worker
    output_func = _output

    for c in kwargs['tgt_data'].keys():
        ref_gts = kwargs['ref_data'][c]['GT']
        tgt_gts = kwargs['tgt_data'][c]['GT']

        ref_mut_num, ref_ind_num, ref_ploidy = ref_gts.shape
        tgt_mut_num, tgt_ind_num, tgt_ploidy = tgt_gts.shape
        ref_hap_num = ref_ind_num*ref_ploidy
        tgt_hap_num = tgt_ind_num*tgt_ploidy
        ref_gts = np.reshape(ref_gts, (ref_mut_num, ref_hap_num))
        tgt_gts = np.reshape(tgt_gts, (tgt_mut_num, tgt_hap_num))

        kwargs['ref_data'][c]['GT'] = ref_gts
        kwargs['tgt_data'][c]['GT'] = tgt_gts
        windows = create_windows(kwargs['tgt_data'][c]['POS'], c, win_step, win_len)


def _process_sstar(win_step, win_len, output, thread, **kwargs):
    """
    Description:
        Processes genotype data with the S* statistic only.

    Arguments:
        win_step int: Length of sliding windows.
        win_len int: Step size of sliding windows.
        output str: Name of the output file.
        thread int: Number of threads.

    Keyword arguments:
        ref_data dict: Dictionary containing data from the reference population.
        tgt_data dict: Dictionary containing data from the target population.
        ref_samples list: List containing sample information from the reference population.
        tgt_samples list: List containing sample information from the target population.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
    """
    header = 'chrom\tstart\tend\tsample\tS*_score\tS*_SNP_number\tS*_SNPs'
    worker_func = _sstar_worker
    output_func = _sstar_output

    for c in kwargs['tgt_data'].keys():
        # Remove variants observed in the reference population
        # Assume 1 is the derived allele
        variants_not_in_ref = np.sum(kwargs['ref_data'][c]['GT'].is_hom_ref(), axis=1) == len(kwargs['ref_samples'])
        kwargs['tgt_data'] = filter_data(kwargs['tgt_data'], c, variants_not_in_ref)
        kwargs['tgt_data'][c]['GT'] = np.sum(kwargs['tgt_data'][c]['GT'], axis=2)
        windows = create_windows(kwargs['tgt_data'][c]['POS'], c, win_step, win_len)


def _output(output, header, samples, res):
    """
    Description:
        Outputs features.

    Arguments:
        output str: Name of the output file.
        header str: Header line for the output file.
        samples list: List containing sample information.
        res list: List containing results.
    """
    with open(output, 'w') as o:
        o.write(header+'\n')
        for item in res:
            chr_name = item[0]
            start = item[1]
            end = item[2]
            spectra = item[3]
            min_ref_dists = item[4]
            tgt_dists = item[5]
            mean_tgt_dists = item[6]
            var_tgt_dists = item[7]
            skew_tgt_dists = item[8]
            kurtosis_tgt_dists = item[9]
            pvt_mut_nums = item[10]
            sstar_scores = item[11]
            for i in range(len(samples)*2):
                ind_name = samples[int(i/2)]
                hap_name = f'hap_{i%2}'
                spectrum = "\t".join([str(x) for x in spectra[i]])
                min_ref_dist = min_ref_dists[i]
                tgt_dist = "\t".join([str(x) for x in tgt_dists[i]])
                mean_tgt_dist = mean_tgt_dists[i]
                var_tgt_dist = var_tgt_dists[i]
                skew_tgt_dist = skew_tgt_dists[i]
                kurtosis_tgt_dist = kurtosis_tgt_dists[i]
                pvt_mut_num = pvt_mut_nums[i]
                sstar_score = sstar_scores[i]
                o.write(f'{chr_name}\t{start}\t{end}\t{ind_name}\t{hap_name}\t{spectrum}\t{tgt_dist}\t{mean_tgt_dist}\t{var_tgt_dist}\t{skew_tgt_dist}\t{kurtosis_tgt_dist}\t{min_ref_dist}\t{sstar_score}\t{pvt_mut_num}\n')


def _sstar_output(output, header, samples, res):
    """
    Description:
        Outputs results from sstar.

    Arguments:
        output str: Name of the output file.
        header str: Header line for the output file.
        samples list: List containing sample information.
        res list: List containing results.
    """
    with open(output, 'w') as o:
        o.write(header+'\n')
        for item in res:
            chr_name = item[0]
            start = item[1]
            end = item[2]
            sstar_scores = item[3]
            sstar_snp_nums = item[4]
            haplotypes = item[5]
            for i in range(len(samples)):
                ind_name = samples[i]
                o.write(f'{chr_name}\t{start}\t{end}\t{ind_name}\t{sstar_scores[i]}\t{sstar_snp_nums[i]}\t{haplotypes[i]}\n')


if __name__ == '__main__':
    process_data(vcf_file="examples/data/real_data/sstar.example.biallelic.snps.vcf.gz", ref_ind_file="examples/data/ind_list/ref.ind.list", tgt_ind_file="examples/data/ind_list/tgt.ind.list", anc_allele_file=None, feature_file="examples/pre-trained/test.features.yaml", output="sstar/test.preprocess.out", win_len=50000, win_step=10000, thread=1)
