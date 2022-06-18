# Apache License Version 2.0
# Copyright 2022 Xin Huang
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

import allel
import gzip
import os
import numpy as np
import pandas as pd
from multiprocessing import Process, Queue
from sstar.utils import read_data, py2round, read_mapped_region_file, cal_matchpct

#@profile
def cal_match_pct(vcf, ref_ind_file, tgt_ind_file, src_ind_file, anc_allele_file, output, thread, score_file, mapped_region_file):
    """
    Description:
        Calculate p-values for S* haplotypes in the target population with source genomes.

    Arguments:
        vcf str: Name of the VCF file containing genotypes.
        src_vcf str: Name of the VCF file containing genotypes from source populations.
        ref_ind_file str: Name of the file containing sample information from reference populations.
        tgt_ind_file str: Name of the file containing sample information from target populations.
        src_ind_file str: Name of the file containing sample information from source populations.
        anc_allele_file str: Name of the file containing ancestral allele information.
        output str: Name of the output file.
        thread int: Number of threads.
        score_file str: Name of the file containing S* scores calculated by `s-star score`.
        mapped_region_file str: Name of the BED file containing mapped regions.
    """

    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(vcf, ref_ind_file, tgt_ind_file, src_ind_file, anc_allele_file)

    res = []
    chr_names = ref_data.keys()

    mapped_intervals = read_mapped_region_file(mapped_region_file)
    data, windows, samples = _read_score_file(score_file, chr_names, tgt_samples)
    sample_size = len(samples)

    header = 'chrom\tstart\tend\tsample\tmatch_rate\tsrc_sample'

    if thread > 1: thread = min(os.cpu_count()-1, sample_size, thread)    
    res = _cal_tgt_match_pct_manager(data, mapped_intervals, samples, tgt_samples, src_samples, tgt_data, src_data, sample_size, thread)

    with open(output, 'w') as o:
        o.write(header+"\n")
        o.write("\n".join(res)+"\n")

#@profile
def _read_score_file(score_file, chr_names, tgt_samples):
    """
    Description:
        Helper function for reading the file generated by `sstar score`.

    Arguments:
        score_file str: Name of the file containing S* scores generated by `sstar score`.
        chr_names list: List containing names of chromosomes for analysis.
        tgt_samples list: List containing names of samples from the target population for analysis.

    Returns:
        data dict: Dictionary containing S* for analysis.
        windows dict: Dictionary containing windows for analysis.
        header str: Header from the file generated by `sstar score`.
        samples list: List containing names of samples in the target population for analysis.
    """

    data = dict()
    windows = dict()
    for c in chr_names:
        windows[c] = []
    samples = []
    with open(score_file, 'r') as f:
        header = f.readline().rstrip()
        for line in f.readlines():
            line = line.rstrip()
            elements = line.split("\t")
            chr_name = elements[0]
            win_start = elements[1]
            win_end = elements[2]
            sample = elements[3]
            if sample not in tgt_samples: continue
            if elements[6] == 'NA': continue
            if sample not in data.keys(): 
                data[sample] = []
                samples.append(sample)
            data[sample].append(line)
            windows[c].append((int(win_start), int(win_end)))

    return data, windows, samples

def _cal_tgt_match_pct_manager(data, mapped_intervals, samples, tgt_samples, src_samples, tgt_data, src_data, sample_size, thread):
    """
    Description:
        Manager function to calculate match percents in target populations using multiprocessing.

    Arguments:
        data dict: Lines from the output file created by `sstar score`.
        mapped_intervals dict: Dictionary of tuples containing mapped regions across the genome.
        sample list: Sample information for individuals needed to be estimated match percents.
        tgt_samples list: Sample information from target populations.
        src_samples list: Sample information from source populations.
        tgt_data dict: Genotype data from target populations.
        src_data dict: Genotype data from source populations.
        sample_size int: Number of individuals analyzed.
        thread int: Number of threads.

    Returns:
        res list: Match percents for target populations.
    """

    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    res = []
    in_queue, out_queue = Queue(), Queue()
    workers = [Process(target=_cal_tgt_match_pct_worker, args=(in_queue, out_queue, mapped_intervals, tgt_data, src_data, src_samples, len(tgt_samples))) for ii in range(thread)]

    for t in samples:
        index = tgt_samples.index(t)
        in_queue.put((index, data[t]))
    
    try:
        for worker in workers:
            worker.start()
        for s in range(sample_size):
            item = out_queue.get()
            if item != '': res.append(item)
        for worker in workers:
            worker.terminate()
    finally:
        for worker in workers:
            worker.join()

    return res

def _cal_tgt_match_pct_worker(in_queue, out_queue, mapped_intervals, tgt_data, src_data, src_samples, sample_size):
    """
    Description:
        Worker function to calculate match percents in target populations.

    Arguments:
        in_queue multiprocessing.Queue: multiprocessing.Queue instance to receive parameters from the manager.
        out_queue multiprocessing.Queue: multiprocessing.Queue instance to send results back to the manager.
        mapped_intervals dict: Dictionary of tuples containing mapped regions across the genome.
        tgt_data dict: Genotype data from target populations.
        src_data dict: Genotype data from source populations.
        src_samples list: List containing sample information for source populations.
        sample_size int: Number of individuals analyzed.
    """

    while True:
        index, data = in_queue.get()
        res = _cal_match_pct_ind(data, index, mapped_intervals, tgt_data, src_data, src_samples, sample_size )
        out_queue.put("\n".join(res))

#@profile
def _cal_match_pct_ind(data, tgt_ind_index, mapped_intervals, tgt_data, src_data, src_samples, sample_size):
    """
    Description:
        Helper function for calculating p-values in individuals.

    Arguments:
        data dict: Dictionary containing S* for analysis.
        tgt_ind_index int: Index of the target individual for analysis.
        mapped_intervals dict: Dictionary of tuples containing mapped regions across the genome.
        tgt_data dict: Genotype data from target populations.
        src_data dict: Genotype data from source populations.
        src_samples list: List of samples from source populations.
        sample_size int: Number of individuals analyzed.

    Returns:
        res list: List containing estimated p-values and other statistics.
    """
    res = []
    for line in data:
        elements = line.split("\t")
        chr_name = elements[0]
        win_start, win_end = elements[1], elements[2]
        sample = elements[3]

        s_star_snps = elements[-1].split(",")
        s_start, s_end = s_star_snps[0], s_star_snps[-1]
        key1 = win_start+'-'+win_end
        key2 = s_start+'-'+s_end

        for src_ind_index in range(len(src_samples)):
            src_sample = src_samples[src_ind_index]
            hap1_res = cal_matchpct(chr_name, mapped_intervals, tgt_data, src_data, tgt_ind_index, src_ind_index, 0, int(win_start), int(win_end), sample_size)
            hap2_res = cal_matchpct(chr_name, mapped_intervals, tgt_data, src_data, tgt_ind_index, src_ind_index, 1, int(win_start), int(win_end), sample_size)
            
            hap1_match_pct = hap1_res[-1]
            hap2_match_pct = hap2_res[-1]
            hap_match_pct = 'NA'

            if (hap1_match_pct != 'NA') and (hap2_match_pct != 'NA'): hap_match_pct = (hap1_match_pct + hap2_match_pct) / 2

            line = f'{chr_name}\t{win_start}\t{win_end}\t{sample}\t{hap_match_pct}\t{src_sample}'
            res.append(line)
    return res