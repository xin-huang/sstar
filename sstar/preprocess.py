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
    windows = {}
    for c in chr_names:
        # Remove variants observed in the reference population
        # Assume 1 is the derived allele
        variants_not_in_ref = np.sum(ref_data[c]['GT'].is_hom_ref(), axis=1) == len(ref_samples)
        tgt_data = filter_data(tgt_data, c, variants_not_in_ref)
        if c not in windows.keys(): windows[c] = []
        else: windows[c] = create_windows(tgt_data[c]['POS'], win_step, win_len)

    if process_archie: _process_archie()
    else: _process_sstar()


def _process_archie():
    """
    """
    pass


def _process_sstar():
    """
    """
    pass


def _manager(worker_func, header):
    """
    """
    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()


def _archie_worker():
    """
    """
    pass


def _sstar_worker():
    """
    """
    pass
