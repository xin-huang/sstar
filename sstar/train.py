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


import pandas as pd
from multiprocessing import Process, Queue
from sstar.preprocess import process_data
from sstar.models import LogisticRegression, ExtraTrees, Sstar
from sstar.utils import multiprocessing_manager


def train(nrep, seq_len, thread, training_data_prefix, training_data_dir, model_file,
          feature_config, is_phased, ploidy, output_genotypes, archaic_prop, not_archaic_prop, algorithm=None):
    """
    """
    multiprocessing_manager(worker_func=_preprocess_worker, nrep=nrep, thread=thread,
                            training_data_prefix=training_data_prefix, training_data_dir=training_data_dir,
                            win_len=seq_len, win_step=seq_len, feature_config=feature_config,
                            is_phased=is_phased, ploidy=ploidy, output_genotypes=output_genotypes,
                            seq_len=seq_len, archaic_prop=archaic_prop, not_archaic_prop=not_archaic_prop)

    feature_df = pd.DataFrame()
    for i in range(nrep):
        feature_file = training_data_dir + '/' + str(i) + '/' + training_data_prefix + f'.{i}.features'
        df = pd.read_csv(feature_file, sep="\t")
        feature_df = pd.concat([feature_df, df])

    all_feature_file = training_data_dir + '/' + training_data_prefix + '.training.all.features'

    feature_df.to_csv(all_feature_file, sep="\t", index=False)
    feature_df = feature_df.drop(columns=['chrom', 'start', 'end', 'sample'])

    if algorithm == 'logistic_regression':
        model = LogisticRegression()
    elif algorithm == 'extra_trees':
        model = ExtraTrees()
    elif (algorithm == 'sstar') or (algorithm is None):
        model = Sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')

    model.train(feature_df, model_file)


def _preprocess_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        rep = in_queue.get()

        vcf_file = kwargs['training_data_dir'] + '/' + str(rep) + '/' + kwargs['training_data_prefix'] + f'.{rep}.vcf'
        bed_file = kwargs['training_data_dir'] + '/' + str(rep) + '/' + kwargs['training_data_prefix'] + f'.{rep}.true.tracts.bed'
        ref_ind_file = kwargs['training_data_dir'] + '/' + str(rep) + '/' + kwargs['training_data_prefix'] + f'.{rep}.ref.ind.list'
        tgt_ind_file = kwargs['training_data_dir'] + '/' + str(rep) + '/' + kwargs['training_data_prefix'] + f'.{rep}.tgt.ind.list'
        feature_file = kwargs['training_data_dir'] + '/' + str(rep) + '/' + kwargs['training_data_prefix'] + f'.{rep}.features'
        label_file = kwargs['training_data_dir'] + '/' + str(rep) + '/' + kwargs['training_data_prefix'] + f'.{rep}.true.tracts.label'

        process_data(vcf_file=vcf_file, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file,
                     anc_allele_file=None, feature_config=kwargs['feature_config'], thread=1,
                     is_phased=kwargs['is_phased'], ploidy=kwargs['ploidy'], output_genotypes=kwargs['output_genotypes'],
                     output_dir=kwargs['training_data_dir']+'/'+str(rep), output_prefix=kwargs['training_data_prefix']+f'.{rep}',
                     win_len=kwargs['win_len'], win_step=kwargs['win_step'])

        feature_df = pd.read_csv(feature_file, sep="\t")

        try: 
            true_tract_df = pd.read_csv(bed_file, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            feature_df['label'] = 0.0
        else:
            true_tract_df.columns = ['chrom', 'start', 'end', 'sample']
            true_tract_df['len'] = true_tract_df['end'] - true_tract_df['start']
            true_tract_df = true_tract_df.groupby(by=['sample'])['len'].sum().reset_index()
            true_tract_df['prop'] = true_tract_df['len'] / kwargs['seq_len']
            true_tract_df['label'] = true_tract_df.apply(lambda row: _add_label(row, kwargs['archaic_prop'], kwargs['not_archaic_prop']), axis=1)
            true_tract_df.to_csv(label_file, sep="\t", index=False)
            feature_df = feature_df.merge(true_tract_df.drop(columns=['len', 'prop']), 
                                          left_on=['sample'], right_on=['sample'], how='left').fillna(0)
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
    train(nrep=1000, thread=2, training_data_prefix='test', training_data_dir='./sstar/test', 
          seq_len=50000, feature_config='examples/pre-trained/test.features.yaml', model_file="./sstar/test/test.lr.model",
          is_phased=True, ploidy=2, output_genotypes=False,
          archaic_prop=0.7, not_archaic_prop=0.3, algorithm='logistic_regression')
