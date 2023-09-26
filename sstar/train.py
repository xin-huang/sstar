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
    # simulate data
    multiprocessing_manager(worker_func=_simulation_worker, nrep=nrep, thread=thread,
             demo_model_file=demo_model_file, nref=nref, ntgt=ntgt, 
             ref_id=ref_id, tgt_id=tgt_id, src_id=src_id, 
             seq_len=seq_len, mut_rate=mut_rate, rec_rate=rec_rate, 
             output_prefix=output_prefix, output_dir=output_dir, seed=seed)

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


def _simulation_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        rep = in_queue.get()

        demo_graph = demes.load(kwargs['demo_model_file'])
        demography = msprime.Demography.from_demes(demo_graph)
        samples = [
            msprime.SampleSet(kwargs['nref'], ploidy=2, population=kwargs['ref_id']),
            msprime.SampleSet(kwargs['ntgt'], ploidy=2, population=kwargs['tgt_id']),
        ]

        ts = msprime.sim_ancestry(
            recombination_rate=kwargs['rec_rate'],
            sequence_length=kwargs['seq_len'],
            samples=samples,
            demography=demography,
            record_migrations=True,
            random_seed=kwargs['seed'],
        )
        ts = msprime.sim_mutations(ts, rate=kwargs['mut_rate'], random_seed=kwargs['seed'], model=msprime.BinaryMutationModel())
        true_tracts = _get_true_tracts(ts, kwargs['tgt_id'], kwargs['src_id'])

        os.makedirs(os.path.join(kwargs['output_dir'], str(rep)), exist_ok=True)
        ts_file  = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.ts'
        vcf_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.vcf'
        bed_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.true.tracts.bed'
        ref_ind_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.ref.ind.list'
        tgt_ind_file = kwargs['output_dir'] + '/' + str(rep) + '/' + kwargs['output_prefix'] + f'.{rep}.tgt.ind.list'

        _create_ref_tgt_file(kwargs['nref'], kwargs['ntgt'], ref_ind_file, tgt_ind_file)

        ts.dump(ts_file)
        with open(vcf_file, 'w') as o:
            ts.write_vcf(o)
       
        df = pd.DataFrame()
        for n in sorted(true_tracts.keys()):
            true_tracts[n].sort(key=lambda x:(x[0], x[1], x[2]))
            df2 = pd.DataFrame(true_tracts[n], columns=['chr', 'start', 'end', 'hap', 'ind'])
            df = pd.concat([df, df2])

        df.drop_duplicates(keep='first').to_csv(bed_file, sep="\t", header=False, index=False)

        out_queue.put(rep)


def _create_ref_tgt_file(nref, ntgt, ref_ind_file, tgt_ind_file, identifier="tsk_"):
    """
    Description:
        Helper function that creates approriate reference and target individual files.

    Arguments:
        nref int: number of reference individuals
        ntgt int: number of target individuals
        ref_ind_file str: Name of the reference individuals output file
        tgt_ind_file str: Name of the target individuals output file
        identifier str: string to prepend at the beginning of the individual number
    """
    with open(ref_ind_file, 'w') as f:
        for i in range(nref):
            f.write(identifier + str(i) + "\n")

    with open(tgt_ind_file, 'w') as f:
        for i in range(i+1, nref + ntgt):
            f.write(identifier + str(i) + "\n")


def _get_true_tracts(ts, tgt_id, src_id):
    """
    Description:
        Helper function to obtain ground truth introgressed tracts from tree-sequence.

    Arguments:
        ts tskit.TreeSqueuece: Tree-sequence containing ground truth introgressed tracts.
        tgt_id str: Name of the target population. 
        src_id str: Name of the source population.
    """
    tracts = {}
    introgression = []

    for p in ts.populations():
        source_id = [p.id for p in ts.populations() if p.metadata['name']==src_id][0]
        target_id = [p.id for p in ts.populations() if p.metadata['name']==tgt_id][0]

    for i in range(ts.num_samples):
        node = ts.node(i)
        if node.population == target_id: tracts[node.id] = []

    for m in ts.migrations():
        if m.dest == source_id: introgression.append(m)

    for i in introgression:
        for t in ts.trees():
            if i.left > t.interval[0]: continue
            if i.right <= t.interval[0]: break
            for n in tracts.keys():
                if t.is_descendant(n, i.node): tracts[n].append([1, int(i.left), int(i.right), f'hap_{int(n%2)}', f'tsk_{ts.node(n).individual}'])

    return tracts


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
