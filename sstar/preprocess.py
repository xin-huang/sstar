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


def process_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, feature_file, output, win_len, win_step, thread, phased):
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
        phased bool: If True, use phased genotypes; otherwise, use unphased genotypes.
    """
    with open(feature_file, 'r') as f:
        features = yaml.safe_load(f)

    print(features)

    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(vcf_file, ref_ind_file, tgt_ind_file, None, anc_allele_file, phased)

    for c in tgt_data.keys():
        windows = create_windows(tgt_data[c]['POS'], c, win_step, win_len)

    res = multiprocessing_manager(worker_func=preprocess_worker, nrep=len(windows), thread=thread, windows=windows, ref_data=ref_data, tgt_data=tgt_data, features=features)

    # x[0]: the chromosome name in number
    # x[1]: the start of the window
    # x[2]: the end of the window
    res.sort(key=lambda x: (x[0], x[1], x[2]))

    print(res)

    #if 'genotypes' in features['features']: 
    #    for c in chr_names:
    #        print(tgt_data[c]['GT'][0])

    #if only_sstar: 
    #    _process_sstar(win_step, win_len, output, thread, 
    #                   ref_data=ref_data, tgt_data=tgt_data, ref_samples=ref_samples, tgt_samples=tgt_samples, 
    #                   match_bonus=match_bonus, max_mismatch=max_mismatch, mismatch_penalty=mismatch_penalty)
    #else:
    #    _process(win_step, win_len, output, thread, 
    #             ref_data=ref_data, tgt_data=tgt_data, samples=tgt_samples, 
    #             match_bonus=match_bonus, max_mismatch=max_mismatch, mismatch_penalty=mismatch_penalty)


def preprocess_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        chr_name, start, end, ref_gts, tgt_gts, pos = in_queue.get()

        if 'genotypes' in kwargs['features']['features']:
            out_queue.put((chr_name, start, end, ref_gts, tgt_gts))
        #spectra = cal_n_ton(tgt_gts)
        #min_ref_dists = cal_ref_dist(ref_gts, tgt_gts)
        #tgt_dists, mean_tgt_dists, var_tgt_dists, skew_tgt_dists, kurtosis_tgt_dists = cal_tgt_dist(tgt_gts)

        #variants_not_in_ref = np.sum(ref_gts, axis=1)==0
        #sub_ref_gts = ref_gts[variants_not_in_ref]
        #sub_tgt_gts = tgt_gts[variants_not_in_ref]
        #sub_pos = pos[variants_not_in_ref]

        #pvt_mut_nums = cal_pvt_mut_num(sub_ref_gts, sub_tgt_gts)
        #sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 'archie', match_bonus, max_mismatch, mismatch_penalty)
        #out_queue.put((chr_name, start, end, 
        #               spectra, min_ref_dists, tgt_dists, 
        #               mean_tgt_dists, var_tgt_dists, skew_tgt_dists, 
        #               kurtosis_tgt_dists, pvt_mut_nums, sstar_scores))


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

    _manager(windows, output, thread, header, worker_func, output_func, 
             ref_data=kwargs['ref_data'], tgt_data=kwargs['tgt_data'], samples=kwargs['samples'],
             match_bonus=kwargs['match_bonus'], max_mismatch=kwargs['max_mismatch'], mismatch_penalty=kwargs['mismatch_penalty'])


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

    _manager(windows, output, thread, header, worker_func, output_func,
             ref_data=None, tgt_data=kwargs['tgt_data'], samples=kwargs['tgt_samples'], 
             match_bonus=kwargs['match_bonus'], max_mismatch=kwargs['max_mismatch'], mismatch_penalty=kwargs['mismatch_penalty'])


def _manager(windows, output, thread, header, worker_func, output_func, **kwargs):
    """
    Description:
        Manages multi-processing jobs.

    Arguments:
        windows list: List containing sliding windows across the genome.
        output str: Name of the output file.
        thread int: Number of threads.
        header str: Header line for the output file.
        worker_func func: Worker function.
        output_func func: Output function.

    Keyword arguments:
        ref_data dict: Dictionary containing data from the reference population.
        tgt_data dict: Dictionary containing data from the target population.
        samples list: List containing sample information from the target population.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
    """
    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    res = []
    in_queue, out_queue = Queue(), Queue()
    workers = [Process(target=worker_func, args=(in_queue, out_queue, kwargs['match_bonus'], kwargs['max_mismatch'], kwargs['mismatch_penalty'])) for i in range(thread)]

    for i in range(len(windows)):
        chr_name, start, end = windows[i]
        tgt_gts = kwargs['tgt_data'][chr_name]['GT']
        pos = kwargs['tgt_data'][chr_name]['POS']
        idx = (pos>start)*(pos<=end)
        sub_ref_gts = None
        sub_tgt_gts = tgt_gts[idx]
        if kwargs['ref_data'] is not None:
            ref_gts = kwargs['ref_data'][chr_name]['GT']
            sub_ref_gts = ref_gts[idx]
        sub_pos = pos[idx]
        in_queue.put((chr_name, start, end, sub_ref_gts, sub_tgt_gts, sub_pos))

    try:
        for w in workers: w.start()

        for i in range(len(windows)):
            item = out_queue.get()
            if item != '': res.append(item)

        for w in workers: w.terminate()
    finally:
        for w in workers: w.join()

    res.sort(key=lambda x: (x[0], x[1], x[2]))
    output_func(output, header, kwargs['samples'], res)


def _worker(in_queue, out_queue, match_bonus, max_mismatch, mismatch_penalty):
    """
    Description:
        Calculates statistics.

    Arguments:
        in_queue multiprocessing.Queue: multiprocessing.Queue instance to receive parameters from the manager.
        out_queue multiprocessing.Queue: multiprocessing.Queue instance to send results to the manager.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
    """
    while True:
        chr_name, start, end, ref_gts, tgt_gts, pos = in_queue.get()
        spectra = cal_n_ton(tgt_gts)
        min_ref_dists = cal_ref_dist(ref_gts, tgt_gts)
        tgt_dists, mean_tgt_dists, var_tgt_dists, skew_tgt_dists, kurtosis_tgt_dists = cal_tgt_dist(tgt_gts)

        variants_not_in_ref = np.sum(ref_gts, axis=1)==0
        sub_ref_gts = ref_gts[variants_not_in_ref]
        sub_tgt_gts = tgt_gts[variants_not_in_ref]
        sub_pos = pos[variants_not_in_ref]

        pvt_mut_nums = cal_pvt_mut_num(sub_ref_gts, sub_tgt_gts)
        sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 'archie', match_bonus, max_mismatch, mismatch_penalty)
        out_queue.put((chr_name, start, end, 
                       spectra, min_ref_dists, tgt_dists, 
                       mean_tgt_dists, var_tgt_dists, skew_tgt_dists, 
                       kurtosis_tgt_dists, pvt_mut_nums, sstar_scores))


def _sstar_worker(in_queue, out_queue, match_bonus, max_mismatch, mismatch_penalty):
    """
    Description:
        Calculates S* scores for sstar.

    Arguments:
        in_queue multiprocessing.Queue: multiprocessing.Queue instance to receive parameters from the manager.
        out_queue multiprocessing.Queue: multiprocessing.Queue instance to send results to the manager.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
    """
    while True:
        chr_name, start, end, ref_gts, tgt_gts, pos = in_queue.get()
        sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(tgt_gts, pos, 'vernot2016', match_bonus, max_mismatch, mismatch_penalty)
        out_queue.put((chr_name, start, end, sstar_scores, sstar_snp_nums, haplotypes))


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
    process_data(vcf_file="examples/data/real_data/sstar.example.biallelic.snps.vcf.gz", ref_ind_file="examples/data/ind_list/ref.ind.list", tgt_ind_file="examples/data/ind_list/tgt.ind.list", anc_allele_file=None, feature_file="examples/pre-trained/test.features.yaml", output="sstar/test.preprocess.out", win_len=50000, win_step=10000, thread=1, phased=False)
