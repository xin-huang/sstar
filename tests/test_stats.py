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


import pytest
import numpy as np
from sstar.stats import *


@pytest.fixture
def data():
    pytest.ref_gts = np.array([
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [1,1,1,1],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [1,1,1,1],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [1,1,1,1],
        [0,0,0,0],
        [0,0,0,0],
        [1,1,1,1],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
    ])

    pytest.tgt_gts = np.array([
        [1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [1,1,1,0,1,1,1,1],
        [0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,1,0,0,0,0,0],
        [1,1,0,1,1,1,1,1],
        [1,0,0,0,0,0,0,0],
        [1,1,1,0,1,1,1,1],
        [1,1,0,1,1,1,1,1],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [1,1,0,1,1,1,1,1],
        [1,1,0,1,1,1,1,1],
        [0,0,0,1,0,0,0,0],
    ])

    pytest.pos = np.array([
       2309,  7879, 11484, 16249, 17324, 19064, 19124, 23559, 25354,
       26654, 29724, 30769, 31319, 37199, 38009, 39444, 40809, 45079,
       48989
    ])

    pytest.spectra = np.array([
       [0, 2, 0, 0, 0, 0, 0, 6, 0],
       [0, 1, 0, 0, 0, 0, 0, 6, 0],
       [0, 2, 0, 0, 0, 0, 0, 2, 0],
       [0, 4, 0, 0, 0, 0, 0, 4, 0],
       [0, 1, 0, 0, 0, 0, 0, 6, 0],
       [0, 1, 0, 0, 0, 0, 0, 6, 0],
       [0, 0, 0, 0, 0, 0, 0, 6, 0],
       [0, 0, 0, 0, 0, 0, 0, 6, 0]
    ])

    pytest.min_ref_dists = np.array([
       2.82842712, 2.64575131, 2.44948974, 3.16227766, 
       2.64575131, 2.64575131, 2.44948974, 2.44948974
    ])

    pytest.tgt_dists = np.array([
        [0.        , 1.41421356, 1.41421356, 1.73205081, 1.73205081,
         1.73205081, 2.82842712, 2.82842712],
        [0.        , 1.        , 1.        , 1.41421356, 1.41421356,
         1.73205081, 2.64575131, 2.64575131],
        [0.        , 2.44948974, 2.44948974, 2.64575131, 2.64575131,
         2.64575131, 2.82842712, 3.46410162],
        [0.        , 2.44948974, 2.44948974, 2.64575131, 2.64575131,
         2.64575131, 2.82842712, 3.46410162],
        [0.        , 1.        , 1.        , 1.41421356, 1.41421356,
         1.73205081, 2.64575131, 2.64575131],
        [0.        , 1.        , 1.        , 1.41421356, 1.41421356,
         1.73205081, 2.64575131, 2.64575131],
        [0.        , 0.        , 1.        , 1.        , 1.        ,
         1.41421356, 2.44948974, 2.44948974],
        [0.        , 0.        , 1.        , 1.        , 1.        ,
         1.41421356, 2.44948974, 2.44948974]
    ])

    pytest.mean_tgt_dists = np.array([
        1.71017922, 1.48149757, 2.39109527, 2.39109527, 
        1.48149757, 1.48149757, 1.16414913, 1.16414913
    ])

    pytest.var_tgt_dists = np.array([
        0.70028702, 0.68016495, 0.90766341, 0.90766341, 
        0.68016495, 0.68016495, 0.7697568 , 0.7697568 
    ])

    pytest.skew_tgt_dists = np.array([
        -0.48140513, -0.06763628, -1.77823182, -1.77823182, 
        -0.06763628, -0.06763628,  0.20248351,  0.20248351
    ])

    pytest.kurtosis_tgt_dists = np.array([
        -0.01859721, -0.67536261,  2.16820852,  2.16820852, 
        -0.67536261, -0.67536261, -1.07216563, -1.07216563
    ])

    pytest.pvt_mut_num = np.array([6, 5, 3, 7, 5, 5, 4, 4])

    pytest.sstar_scores = np.array([51470.0, 86110.0, 36005.0, 33425.0])
    pytest.sstar_snp_num = np.array([6, 10, 6, 4])
    pytest.haplotypes = np.array([
        '2309,25354,26654,29724,40809,45079', 
        '7879,17324,19124,26654,29724,30769,38009,40809,45079,48989', 
        '11484,19064,26654,29724,40809,45079', 
        '26654,29724,40809,45079'
    ])
    pytest.archie_sstar_scores = np.array([25355.0, 25355.0, -np.inf, 25355.0, 25355.0, 25355.0, 25355.0, 25355.0])


def test_cal_n_ton(data):
    spectra = cal_n_ton(pytest.tgt_gts)

    assert np.array_equal(spectra, pytest.spectra)


def test_cal_ref_dist(data):
    min_ref_dists = cal_ref_dist(pytest.ref_gts, pytest.tgt_gts)

    assert np.allclose(min_ref_dists, pytest.min_ref_dists)


def test_cal_tgt_dist(data):
    tgt_dists, mean_tgt_dists, var_tgt_dists, skew_tgt_dists, kurtosis_tgt_dists = cal_tgt_dist(pytest.tgt_gts)

    assert np.allclose(tgt_dists, pytest.tgt_dists)
    assert np.allclose(mean_tgt_dists, pytest.mean_tgt_dists)
    assert np.allclose(var_tgt_dists, pytest.var_tgt_dists)
    assert np.allclose(skew_tgt_dists, pytest.skew_tgt_dists)
    assert np.allclose(kurtosis_tgt_dists, pytest.kurtosis_tgt_dists)


def test_cal_pvt_mut_num(data):
    pvt_mut_num = cal_pvt_mut_num(pytest.ref_gts, pytest.tgt_gts)

    assert np.array_equal(pvt_mut_num, pytest.pvt_mut_num)


def test_cal_sstar(data):
    variants_not_in_ref = np.sum(pytest.ref_gts, axis=1)==0
    sub_ref_gts = pytest.ref_gts[variants_not_in_ref]
    sub_pos = pytest.pos[variants_not_in_ref]

    tgt_gts = np.reshape(pytest.tgt_gts, (pytest.tgt_gts.shape[0],int(pytest.tgt_gts.shape[1]/2),2))
    tgt_gts = np.sum(tgt_gts, axis=2)
    sub_tgt_gts = tgt_gts[variants_not_in_ref]
    sstar_scores, sstar_snp_num, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 'vernot2016', 5000, 5, -10000)

    assert np.array_equal(sstar_scores, pytest.sstar_scores)
    assert np.array_equal(sstar_snp_num, pytest.sstar_snp_num)
    assert np.array_equal(haplotypes, pytest.haplotypes)

    sub_tgt_gts = pytest.tgt_gts[variants_not_in_ref]
    sstar_scores, sstar_snp_num, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 'archie', 5000, 5, -10000)

    assert np.array_equal(sstar_scores, pytest.archie_sstar_scores)
