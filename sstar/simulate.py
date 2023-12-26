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


import os, shutil
import demes, msprime
import pandas as pd
from multiprocessing import Process, Queue
from sstar.utils import multiprocessing_manager
from sstar.preprocess import preprocess


def simulate(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, ploidy,
             feature_config, is_phased, intro_prop, not_intro_prop, rm_sim_data,
             seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed):
    """
    """
    multiprocessing_manager(worker_func=_simulation_worker, nrep=nrep, thread=thread,
                            demo_model_file=demo_model_file, nref=nref, ntgt=ntgt, 
                            ref_id=ref_id, tgt_id=tgt_id, src_id=src_id, ploidy=ploidy,
                            seq_len=seq_len, mut_rate=mut_rate, rec_rate=rec_rate, 
                            feature_config=feature_config, is_phased=is_phased, 
                            intro_prop=intro_prop, not_intro_prop=not_intro_prop,
                            output_prefix=output_prefix, output_dir=output_dir, seed=seed)

    if feature_config is not None:
        feature_df = pd.DataFrame()
        for i in range(nrep):
            df = pd.read_csv(f'{output_dir}/{i}/{output_prefix}.{i}.labeled.features', sep="\t")
            df.insert(0, 'replicate', i)
            feature_df = pd.concat([feature_df, df])

        feature_df.to_csv(f'{output_dir}/{output_prefix}.all.labeled.features', sep="\t", index=False)

    if rm_sim_data is True:
        for i in range(nrep):
            shutil.rmtree(f'{output_dir}/{i}', ignore_errors=True)


def _simulation_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        rep, seed = in_queue.get()

        demo_graph = demes.load(kwargs['demo_model_file'])
        demography = msprime.Demography.from_demes(demo_graph)
        samples = [
            msprime.SampleSet(kwargs['nref'], ploidy=kwargs['ploidy'], population=kwargs['ref_id']),
            msprime.SampleSet(kwargs['ntgt'], ploidy=kwargs['ploidy'], population=kwargs['tgt_id']),
        ]

        ts = msprime.sim_ancestry(
            recombination_rate=kwargs['rec_rate'],
            sequence_length=kwargs['seq_len'],
            samples=samples,
            demography=demography,
            record_migrations=True,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(ts, rate=kwargs['mut_rate'], random_seed=seed, model=msprime.BinaryMutationModel())
        truth_tracts = _get_truth_tracts(ts, kwargs['tgt_id'], kwargs['src_id'], kwargs['ploidy'])

        os.makedirs(f'{kwargs["output_dir"]}/{rep}', exist_ok=True)
        ts_file  = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.ts'
        vcf_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.vcf'
        bed_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.truth.tracts.bed'
        ref_ind_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.ref.ind.list'
        tgt_ind_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.tgt.ind.list'

        _create_ref_tgt_file(kwargs['nref'], kwargs['ntgt'], ref_ind_file, tgt_ind_file)

        ts.dump(ts_file)
        with open(vcf_file, 'w') as o: ts.write_vcf(o)
       
        df = pd.DataFrame()
        for n in sorted(truth_tracts.keys()):
            truth_tracts[n].sort(key=lambda x:(x[0], x[1], x[2]))
            df2 = pd.DataFrame(truth_tracts[n], columns=['chr', 'start', 'end', 'sample'])
            df = pd.concat([df, df2])

        df.drop_duplicates(keep='first').to_csv(bed_file, sep="\t", header=False, index=False)

        if kwargs['feature_config'] is not None:
            preprocess(vcf_file=vcf_file, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file, anc_allele_file=None,
                       feature_config=kwargs["feature_config"], is_phased=kwargs["is_phased"],
                       ploidy=kwargs["ploidy"], output_dir=f'{kwargs["output_dir"]}/{rep}',
                       output_prefix=f'{kwargs["output_prefix"]}.{rep}', win_len=kwargs["seq_len"], win_step=kwargs["seq_len"], thread=1)

            feature_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.features'
            output = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.labeled.features'

            _label(feature_file=feature_file, truth_tract_file=bed_file, output=output,
                   seq_len=kwargs["seq_len"], intro_prop=kwargs["intro_prop"], not_intro_prop=kwargs["not_intro_prop"])

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


def _get_truth_tracts(ts, tgt_id, src_id, ploidy):
    """
    Description:
        Helper function to obtain ground truth introgressed tracts at the haploid level from tree-sequence.

    Arguments:
        ts tskit.TreeSqueuece: Tree-sequence containing ground truth introgressed tracts.
        tgt_id str: Name of the target population. 
        src_id str: Name of the source population.
        ploidy int: Ploidy of the genomes.
    """
    tracts = {}
    introgression = []

    src_id = [p.id for p in ts.populations() if p.metadata['name']==src_id][0]
    tgt_id = [p.id for p in ts.populations() if p.metadata['name']==tgt_id][0]

    for i in ts.samples(tgt_id): tracts[i] = []

    for m in ts.migrations():
        if (m.dest==src_id) and (m.source==tgt_id):
            # For simulations with a long sequence, large sample size, and/or deep generation time
            # This function may become slow
            # May parallelize this function when necessary
            for t in ts.trees():
                # Tree-sequences are sorted by the left ends of the intervals
                # Can skip those tree-sequences are not overlapped with the interval of i.
                if m.left >= t.interval.right: continue
                if m.right <= t.interval.left: break # [l, r)
                for n in tracts.keys():
                    if t.is_descendant(n, m.node):
                        left = m.left if m.left > t.interval.left else t.interval.left
                        right = m.right if m.right < t.interval.right else t.interval.right
                        tracts[n].append([1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}'])

    return tracts


def _label(feature_file, truth_tract_file, seq_len, intro_prop, not_intro_prop, output):
    """
    """
    feature_df = pd.read_csv(feature_file, sep="\t")

    try:
        truth_tract_df = pd.read_csv(truth_tract_file, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        feature_df['label'] = 0.0
    else:
        truth_tract_df.columns = ['chrom', 'start', 'end', 'sample']
        truth_tract_df['len'] = truth_tract_df['end'] - truth_tract_df['start']
        truth_tract_df = truth_tract_df.groupby(by=['sample'])['len'].sum().reset_index()
        truth_tract_df['prop'] = truth_tract_df['len'] / seq_len
        truth_tract_df['label'] = truth_tract_df.apply(lambda row: _add_label(row, intro_prop, not_intro_prop), axis=1)
        feature_df = feature_df.merge(truth_tract_df.drop(columns=['len', 'prop']),
                                      left_on=['sample'], right_on=['sample'], how='left').fillna(0)
    finally:
        feature_df.to_csv(output, sep="\t", index=False)


def _add_label(row, intro_prop, not_intro_prop):
    """
    """
    if row['prop'] > intro_prop: return 1.0
    elif row['prop'] < not_intro_prop: return 0.0
    else: return -1.0


if __name__ == '__main__':
    simulate(demo_model_file="./examples/models/ArchIE_3D19.yaml", nrep=1000, nref=50, ntgt=50, 
             ref_id='Ref', tgt_id='Tgt', src_id='Ghost', ploidy=2, seq_len=50000, mut_rate=1.25e-8, rec_rate=1e-8, thread=2,
             feature_config='./examples/features/archie.features.yaml', is_phased=True, intro_prop=0.7, not_intro_prop=0.3, rm_sim_data=True,
             output_prefix='test', output_dir='./sstar/test', seed=913)
