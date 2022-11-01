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


import os
import numpy as np
from multiprocessing import Process, Queue
from sstar.stats import *
from sstar.utils import read_data, filter_data, create_windows


def process_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, output, win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, process_archie=False):
    """
    """
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(vcf_file, ref_ind_file, tgt_ind_file, None, anc_allele_file)
    chr_names = tgt_data.keys()
    for c in chr_names:
        # Remove variants observed in the reference population
        # Assume 1 is the derived allele
        variants_not_in_ref = np.sum(ref_data[c]['GT'].is_hom_ref(), axis=1) == len(ref_samples)
        tgt_data = filter_data(tgt_data, c, variants_not_in_ref)
        windows = create_windows(tgt_data[c]['POS'], c, win_step, win_len)

    if process_archie: _process_archie()
    else: _process_sstar(windows, output, thread, data=tgt_data, samples=tgt_samples, match_bonus=match_bonus, max_mismatch=max_mismatch, mismatch_penalty=mismatch_penalty)


def _process_archie():
    """
    """
    header = ''
    worker_func = _archie_worker
    output_func = _archie_output
    _manager(windows, output, thread, header, worker_func, output_func, **kwargs)


def _process_sstar(windows, output, thread, **kwargs):
    """
    """
    header = 'chrom\tstart\tend\tsample\tS*_score\tS*_SNP_number\tS*_SNPs'
    worker_func = _sstar_worker
    output_func = _sstar_output
    for c in kwargs['data'].keys():
        kwargs['data'][c]['GT'] = np.sum(kwargs['data'][c]['GT'], axis=2)

    _manager(windows, output, thread, header, worker_func, output_func,
             data=kwargs['data'], samples=kwargs['samples'], 
             match_bonus=kwargs['match_bonus'], max_mismatch=kwargs['max_mismatch'], mismatch_penalty=kwargs['mismatch_penalty'])


def _manager(windows, output, thread, header, worker_func, output_func, **kwargs):
    """
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
        gts = kwargs['data'][chr_name]['GT']
        pos = kwargs['data'][chr_name]['POS']
        idx = (pos>start)*(pos<=end)
        sub_gts = gts[idx]
        sub_pos = pos[idx]
        in_queue.put((chr_name, start, end, sub_gts, sub_pos))

    try:
        for w in workers: w.start()

        for i in range(len(windows)):
            item = out_queue.get()
            if item != '': res.append(item)

        for w in workers: w.terminate()
    finally:
        for w in workers: w.join()

    output_func(output, header, kwargs['samples'], res)


def _archie_worker(in_queue, out_queue, match_bonus, max_mismatch, mismatch_penalty):
    """
    """
    pass


def _sstar_worker(in_queue, out_queue, match_bonus, max_mismatch, mismatch_penalty):
    """
    """
    while True:
        chr_name, start, end, gts, pos = in_queue.get()
        sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(gts, pos, match_bonus, max_mismatch, mismatch_penalty)
        out_queue.put((chr_name, start, end, sstar_scores, sstar_snp_nums, haplotypes))


def _archie_output(output, header, samples, res):
    pass


def _sstar_output(output, header, samples, res):
    with open(output, 'w') as o:
        o.write(header+'\n')
        for item in res:
            chr_name = item[0]
            start = item[1]
            end = item[2]
            sstar_scores = item[3]
            sstar_snp_nums = item[4]
            haplotypes = item[5]
            for i in range(len(sstar_scores)):
                ind_name = samples[i]
                o.write(f'{chr_name}\t{start}\t{end}\t{ind_name}\t{sstar_scores[i]}\t{sstar_snp_nums[i]}\t{haplotypes[i]}\n')


if __name__ == '__main__':
    process_data('../tests/data/test.score.data.vcf', '../tests/data/test.ref.ind.list', '../tests/data/test.tgt.ind.list', None, 'test.out', 50000, 10000, 1, 5000, 5, -10000, process_archie=False)
