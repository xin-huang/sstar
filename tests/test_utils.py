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
import pytest
import numpy as np
from sstar.utils import *
from sstar.utils import _cal_mapped_len, _cal_hap_stats


@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.score.data.vcf"
    pytest.anc_allele = "./tests/data/test.anc.allele.bed"
    pytest.emp_ind_list = "./tests/data/test.empty.ind.list"
    pytest.emp_anc_allele = "./tests/data/test.empty.anc.allele.bed"
    pytest.mapped_regions = "tests/data/test.mapped.region.bed"


def test_parse_inds_file(data):
    ref_ind = parse_ind_file(pytest.ref_ind_list)
    tgt_ind = parse_ind_file(pytest.tgt_ind_list)

    exp_ref_ind = ['ind5', 'ind6']
    exp_tgt_ind = ['ind1', 'ind2', 'ind3', 'ind4']

    assert ref_ind == exp_ref_ind
    assert tgt_ind == exp_tgt_ind

    with pytest.raises(Exception) as e_info:
        emp_ind = parse_ind_file(pytest.emp_ind_list)


def test_read_geno_data(data):
    ref_ind = parse_ind_file(pytest.ref_ind_list)
    d = read_geno_data(pytest.vcf, ref_ind, None, filter_missing=False)

    vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ref_ind)

    assert np.array_equal(ref_ind, vcf['samples'])
    assert np.array_equal(d['21']['POS'], vcf['variants/POS'])
    assert np.array_equal(d['21']['REF'], vcf['variants/REF'])
    assert np.array_equal(d['21']['ALT'], vcf['variants/ALT'])
    assert np.array_equal(d['21']['GT'], vcf['calldata/GT'])


def test_read_data(data):
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(pytest.vcf, pytest.ref_ind_list, pytest.tgt_ind_list, None, None)

    rs = parse_ind_file(pytest.ref_ind_list)
    ts = parse_ind_file(pytest.tgt_ind_list)
    
    assert np.array_equal(rs, ref_samples)
    assert np.array_equal(ts, tgt_samples)

    ref_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=rs)
    tgt_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ts)

    assert np.array_equal(rs, ref_vcf['samples'])
    assert np.array_equal(ts, tgt_vcf['samples'])
    assert np.array_equal(ref_data['21']['POS'], ref_vcf['variants/POS'])
    assert np.array_equal(ref_data['21']['REF'], ref_vcf['variants/REF'])
    assert np.array_equal(ref_data['21']['ALT'], ref_vcf['variants/ALT'])
    assert np.array_equal(ref_data['21']['GT'], ref_vcf['calldata/GT'])
    assert np.array_equal(tgt_data['21']['POS'], tgt_vcf['variants/POS'])
    assert np.array_equal(tgt_data['21']['REF'], tgt_vcf['variants/REF'])
    assert np.array_equal(tgt_data['21']['ALT'], tgt_vcf['variants/ALT'])
    assert np.array_equal(tgt_data['21']['GT'], tgt_vcf['calldata/GT'])


def test_read_anc_allele(data):
    anc_allele = read_anc_allele(pytest.anc_allele)

    exp_anc_allele = {
        '21': {
            2309: 'G', 7879: 'A', 11484: '-', 48989: 'C'
        }
    }

    assert anc_allele == exp_anc_allele

    with pytest.raises(Exception) as e_info:
        anc_allele = read_anc_allele(pytest.emp_anc_allele)


def test_get_ref_alt_allele(data):
    ref_ind = parse_ind_file(pytest.ref_ind_list)
    tgt_ind = parse_ind_file(pytest.tgt_ind_list)

    ref_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=ref_ind)
    tgt_vcf = allel.read_vcf(pytest.vcf, alt_number=1, samples=tgt_ind)

    ref_allele, alt_allele = get_ref_alt_allele(tgt_vcf['variants/REF'], tgt_vcf['variants/ALT'], tgt_vcf['variants/POS'])

    exp_ref_allele = {
                         2309: 'G', 7879: 'C', 11484: 'A', 16249: 'A', 17324: 'G',
                         19064: 'G', 19124: 'G', 23559: 'G', 25354: 'G', 26654: 'G',
                         29724: 'G', 30769: 'C', 31319: 'C', 37199: 'C', 38009: 'C',
                         39444: 'C', 40809: 'C', 45079: 'C', 48989: 'C'
                     }
    exp_alt_allele = {
                         2309: 'A', 7879: 'A', 11484: 'C', 16249: 'C', 17324: 'T', 
                         19064: 'T', 19124: 'A', 23559: 'A', 25354: 'T', 26654: 'C', 
                         29724: 'A', 30769: 'T', 31319: 'T', 37199: 'T', 38009: 'T', 
                         39444: 'T', 40809: 'T', 45079: 'T', 48989: 'T'
                     }
   
    assert ref_allele == exp_ref_allele
    assert alt_allele == exp_alt_allele


def test_check_anc_allele(data):
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(pytest.vcf, pytest.ref_ind_list, pytest.tgt_ind_list, None, pytest.anc_allele)
  
    exp_ref_gt = allel.GenotypeArray([[[0,0], [0,0]],
                                      [[1,1], [1,1]],
                                      [[0,0], [0,0]]], dtype='i1')
    exp_tgt_gt = allel.GenotypeArray([[[1,0], [0,0], [0,0], [0,0]],
                                      [[1,1], [1,0], [1,1], [1,1]],
                                      [[0,0], [0,1], [0,0], [0,0]]], dtype='i1')
    exp_tgt_pos = [2309, 7879, 48989]

    assert np.array_equal(ref_data['21']['GT'], exp_ref_gt)
    assert np.array_equal(tgt_data['21']['GT'], exp_tgt_gt)
    assert np.array_equal(tgt_data['21']['POS'], exp_tgt_pos)


def test_create_windows(data):
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(pytest.vcf, pytest.ref_ind_list, pytest.tgt_ind_list, None, pytest.anc_allele)
    windows = create_windows(tgt_data['21']['POS'], '21', 10000, 50000)

    exp_windows = [('21', 0, 50000), ('21', 10000, 60000), ('21', 20000, 70000), ('21', 30000, 80000), ('21', 40000, 90000)]

    assert windows == exp_windows


def test_cal_mapped_len(data):
    mapped_intervals = read_mapped_region_file(pytest.mapped_regions)
    len1 = _cal_mapped_len(mapped_intervals, '21', 0, 50000)
    len2 = _cal_mapped_len(mapped_intervals, '21', 1000, 4000)
    len3 = _cal_mapped_len(mapped_intervals, '21', 20000, 60000)
    len4 = _cal_mapped_len(mapped_intervals, 21, 20000, 60000)
    len5 = _cal_mapped_len(mapped_intervals, '21', 60000, 100000)
    len6 = _cal_mapped_len(mapped_intervals, '21', 60000, 72000)

    assert len1 == 50000
    assert len2 == 3000
    assert len3 == 30000
    assert len4 == 40000
    assert len5 == 18000
    assert len6 == 2000


def test_cal_matchpct(data):
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data("./tests/data/test.match.rate.data.vcf", "./examples/data/ind_list/ref.ind.list", "./examples/data/ind_list/tgt.ind.list", "./examples/data/ind_list/nean.ind.list", None)

    hap1_match_pct = cal_matchpct('21', None, tgt_data, src_data, 0, 0, 0, 9400000, 9450000, len(tgt_data))[-1]
    hap2_match_pct = cal_matchpct('21', None, tgt_data, src_data, 0, 0, 1, 9400000, 9450000, len(tgt_data))[-1]

    assert hap1_match_pct == 0.083333
    assert hap2_match_pct == 0.068966
