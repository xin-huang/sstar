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

from multiprocessing import Process, Queue
from sstar.utils import parse_ind_file, read_geno_data, read_mapped_region_file, cal_match_pct

def cal_ref_match_pct(vcf, ref_ind_file, src_ind_file, anc_allele_file, output, win_len, win_step, thread, mapped_region_file):
    """
    Description:
        Calculate archaic match percentages from the reference population.

    Arguments:
        vcf str: Name of the VCF file containing genotypes.
        ref_ind_file str: Name of the file containing sample information from reference populations.
        src_ind_file str: Name of the file containing sample information from source populations.
        anc_allele_file str: Name of the file containing ancestral allele information.
        output str: Name of the output file.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        thread int: Number of threads.
        mapped_region_file str: Name of the BED file containing mapped regions.
    """

    ref_samples = parse_ind_file(ref_ind_file)
    src_samples = parse_ind_file(src_ind_file)
    ref_data = read_geno_data(vcf, ref_samples, anc_allele_file, True)
    src_data = read_geno_data(vcf, src_samples, anc_allele_file, False)
    ref_chr_names = ref_data.keys()
    src_chr_names = src_data.keys()
    chr_names = ref_chr_names & src_chr_names
    if len(chr_names) == 0: raise Exception("Cannot find common chromosomes intersected in both "+ref_vcf+" and "+src_vcf)

    mapped_intervals = read_mapped_region_file(mapped_region_file)

    header = "count\thap_len\thap_mapped_len_bin\thap_match_num\thap_tot_num\thap_dSNP_per_site_num\tmatch_pct"

    windows = dict()
    for c in chr_names:
        pos = ref_data[c]['POS']
        windows[c] = []
        win_start = (pos[0] + win_step) // win_step * win_step - win_len
        if win_start < 0: win_start = 0
        last_pos = pos[-1]
        win_end = 0
        while last_pos > win_start:
            win_end = win_start + win_len
            windows[c].append((win_start, win_end))
            win_start += win_step

    res = _cal_ref_match_pct_manager(ref_data, src_data, mapped_intervals, ref_samples, src_samples, windows, thread)

    with open(output, 'w') as o:
        o.write(header+"\n")
        o.write("\n".join(res)+"\n")

def _cal_ref_match_pct_manager(ref_data, src_data, mapped_intervals, ref_samples, src_samples, windows, thread):
    """
    Description:
        Manager function to calculate match percents in reference populations using multiprocessing.

    Arguments:
        ref_data dict: Genotype data from reference populations.
        src_data dict: Genotype data from source populations.
        mapped_intervals dict: Dictionary of tuples containing mapped regions across the genome. 
        ref_samples list: List of samples from reference populations.
        src_samples list: List of samples from source populations.
        windows list: Regions for calculating match percents.
        thread int: Number of threads.

    Returns:
        ref_match_pct dict: Match percents for estimating p-values within a local window.
        res str: Match percents for estimating p-values using whole-genome data.
    """

    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    in_queue, out_queue = Queue(), Queue()
    res = []
    sample_size = len(ref_samples)

    workers = [Process(target=_cal_ref_match_pct_worker, args=(in_queue, out_queue, mapped_intervals, ref_data, src_data, src_samples, sample_size)) for ii in range(thread)]

    for s in range(sample_size):
        in_queue.put((s, windows))

    try:
        for worker in workers:
            worker.start()
        for s in range(sample_size):
            item = out_queue.get()
            if item != '':
                for k in item.keys():
                    count = item[k]
                    res.append(str(count)+"\t"+k)
        for worker in workers:
            worker.terminate()
    finally:
        for worker in workers:
            worker.join()

    return res

def _cal_ref_match_pct_worker(in_queue, out_queue, mapped_intervals, ref_data, src_data, src_samples, sample_size):
    """
    Description:
        Worker function to calculate match percents in reference populations.

    Arguments:
        in_queue multiprocessing.Queue: multiprocessing.Queue instance to receive parameters from the manager.
        out_queue multiprocessing.Queue: multiprocessing.Queue instance to send results back to the manager.
        mapped_intervals dict: Dictionary of tuples containing mapped regions across the genome.
        ref_data dict: Genotype data from reference populations.
        src_data dict: Genotype data from source populations.
        src_samples list: List of samples from source populations.
        sample_size int: Number of individuals analyzed.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
    """

    while True:
        ref_ind_index, windows = in_queue.get()
        ref_match_pct = dict()
        chr_names = ref_data.keys()
        for c in chr_names:
            for w in windows[c]:
                win_start, win_end = w
                for src_ind_index in range(len(src_samples)):
                    hap1_res = cal_match_pct(c, mapped_intervals, ref_data, src_data, ref_ind_index, src_ind_index, 0, win_start, win_end, sample_size)
                    hap2_res = cal_match_pct(c, mapped_intervals, ref_data, src_data, ref_ind_index, src_ind_index, 1, win_start, win_end, sample_size)
                    hap1_res = hap1_res[1:]
                    hap2_res = hap2_res[1:]
                    key1 = "\t".join([str(r) for r in hap1_res])
                    key2 = "\t".join([str(r) for r in hap2_res])
                    if key1 not in ref_match_pct.keys(): ref_match_pct[key1] = 1
                    else: ref_match_pct[key1] += 1
                    if key2 not in ref_match_pct.keys(): ref_match_pct[key2] = 1
                    else: ref_match_pct[key2] += 1
        out_queue.put(ref_match_pct)
