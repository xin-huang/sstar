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
import demes, msprime
import pandas as pd
from multiprocessing import Process, Queue
from sstar.preprocess import process_data
from sstar.models import train_logistic_regression
from sstar.utils import multiprocessing_manager


def train(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, 
          match_bonus, max_mismatch, mismatch_penalty, archaic_prop, not_archaic_prop, algorithm=None, seed=None):
    """
    """
    if algorithm == 'logistic_regression':
        _train_logistic_regression(nrep=nrep, thread=thread, output_prefix=output_prefix, output_dir=output_dir,
                                   seq_len=seq_len, match_bonus=match_bonus, max_mismatch=max_mismatch, 
                                   mismatch_penalty=mismatch_penalty, archaic_prop=archaic_prop, not_archaic_prop=not_archaic_prop)
    elif algorithm == 'extra_trees':
        _train_extra_trees()
    elif (algorithm == 'sstar') or (algorithm is None):
        _train_sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')


def _train_logistic_regression(nrep, thread, output_prefix, output_dir, seq_len, match_bonus, max_mismatch, mismatch_penalty, archaic_prop, not_archaic_prop):
    """
    """
    multiprocessing_manager(worker_func=_preprocess_logistic_regression_worker, nrep=nrep, thread=thread,
             output_prefix=output_prefix, output_dir=output_dir,
             win_len=seq_len, win_step=seq_len, match_bonus=match_bonus, 
             max_mismatch=max_mismatch, mismatch_penalty=mismatch_penalty,
             seq_len=seq_len, archaic_prop=archaic_prop, not_archaic_prop=not_archaic_prop)

    feature_df = pd.DataFrame()
    for i in range(nrep):
        feature_file = output_dir + '/' + str(i) + '/' + output_prefix + f'.{i}.features'
        df = pd.read_csv(feature_file, sep="\t")
        feature_df = pd.concat([feature_df, df])

    all_feature_file = output_dir + '/' + output_prefix + '.all.features'
    model_file = output_dir + '/' + output_prefix + '.logistic.regression.model'

    feature_df.to_csv(all_feature_file, sep="\t", index=False)
    feature_df = feature_df.drop(columns=['chrom', 'start', 'end', 'sample', 'hap'])
    train_logistic_regression(feature_df, model_file)


def _train_extra_trees():
    pass


def _train_sstar():
    pass


def _preprocess_logistic_regression_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        rep = in_queue.get()

        vcf_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.vcf'
        bed_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.true.tracts.bed'
        ref_ind_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.ref.ind.list'
        tgt_ind_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.tgt.ind.list'
        feature_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.features'
        label_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.true.tracts.label'

        process_data(vcf_file=vcf_file, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file,
                     anc_allele_file=None, output=feature_file, thread=1,
                     win_len=kwargs['win_len'], win_step=kwargs['win_step'],
                     match_bonus=kwargs['match_bonus'], max_mismatch=kwargs['max_mismatch'], 
                     mismatch_penalty=kwargs['mismatch_penalty'])

        feature_df = pd.read_csv(feature_file, sep="\t")

        try: 
            true_tract_df = pd.read_csv(bed_file, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            feature_df['label'] = 0.0
        else:
            true_tract_df.columns = ['chr', 'start', 'end', 'hap', 'sample']
            true_tract_df['len'] = true_tract_df['end'] - true_tract_df['start']
            true_tract_df = true_tract_df.groupby(by=['sample', 'hap'])['len'].sum().reset_index()
            true_tract_df['prop'] = true_tract_df['len'] / kwargs['seq_len']
            true_tract_df['label'] = true_tract_df.apply(lambda row: _add_label(row, kwargs['archaic_prop'], kwargs['not_archaic_prop']), axis=1)
            true_tract_df.to_csv(label_file, sep="\t", index=False)
            feature_df = feature_df.merge(true_tract_df.drop(columns=['len', 'prop']), 
                                          left_on=['sample', 'hap'], right_on=['sample', 'hap'], how='left').fillna(0)
        finally:
            feature_df.to_csv(feature_file, sep="\t", index=False)

        out_queue.put(rep)


def _add_label(row, archaic_prop, not_archaic_prop):
    """
    """
    if row['prop'] > archaic_prop: return 1.0
    elif row['prop'] < not_archaic_prop: return 0.0
    else: return -1.0


if __name__ == '__main__':
    train(demo_model_file="./examples/models/BonoboGhost_4K19.yaml", nrep=1000, nref=50, ntgt=50, 
          ref_id='Western', tgt_id='Bonobo', src_id='Ghost', seq_len=50000, mut_rate=1e-8, rec_rate=1e-8, thread=2, 
          output_prefix='test', output_dir='./sstar/test', match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000, 
          archaic_prop=0.7, not_archaic_prop=0.3, algorithm='logistic_regression', seed=913)
