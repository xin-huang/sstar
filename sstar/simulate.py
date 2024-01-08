# Apache License Version 2.0
# Copyright 2024 Xin Huang
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
import pybedtools
import pandas as pd
from multiprocessing import Process, Queue
from sstar.utils import multiprocessing_manager
from sstar.preprocess import preprocess


def simulate(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, ploidy,
             feature_config, is_phased, intro_prop, not_intro_prop, keep_sim_data,
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
        label_df = pd.DataFrame()
        for i in range(nrep):
            df = pd.read_csv(f'{output_dir}/{i}/{output_prefix}.{i}.labels', sep="\t")
            label_df = pd.concat([label_df, df])
        label_df.to_csv(f'{output_dir}/{output_prefix}.all.labels', sep="\t", index=False)

        feature_df = pd.DataFrame()
        for i in range(nrep):
            df = pd.read_csv(f'{output_dir}/{i}/{output_prefix}.{i}.features', sep="\t")
            feature_df = pd.concat([feature_df, df])
        feature_df.to_csv(f'{output_dir}/{output_prefix}.all.features', sep="\t", index=False)
        feature_df = feature_df.merge(label_df, suffixes=('', '_dup'),
                                      left_on=['sample', 'rep'], right_on=['sample', 'rep'], how='left')
        feature_df.drop([i for i in feature_df.columns if 'dup' in i], axis=1, inplace=True)
        feature_df.to_csv(f'{output_dir}/{output_prefix}.all.labeled.features', sep="\t", index=False)

    if keep_sim_data is not True:
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
        if (kwargs["is_phased"] is not True) and (truth_tracts.empty is not True): 
            truth_tracts['sample'] = truth_tracts.apply(lambda row: "_".join(row['sample'].split("_")[0:2]), axis=1)

        os.makedirs(f'{kwargs["output_dir"]}/{rep}', exist_ok=True)
        ts_file  = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.ts'
        vcf_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.vcf'
        bed_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.truth.tracts.bed'
        ref_ind_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.ref.ind.list'
        tgt_ind_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.tgt.ind.list'
        label_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.labels'
        seed_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.seedmsprime'

        _create_ref_tgt_file(kwargs['nref'], kwargs['ntgt'], ref_ind_file, tgt_ind_file)

        ts.dump(ts_file)
        with open(vcf_file, 'w') as o: ts.write_vcf(o)

        with open(seed_file, 'w') as o: o.write(f'{seed}\n')

        open(bed_file, 'w').close()
        for s in truth_tracts['sample'].unique():
            sample_tracts = truth_tracts[truth_tracts['sample'] == s]
            sample_tracts = pybedtools.BedTool.from_dataframe(sample_tracts).sort().merge().to_dataframe()
            sample_tracts['sample'] = s
            sample_tracts.to_csv(bed_file, sep="\t", mode='a', header=False, index=False)

        if kwargs['feature_config'] is not None:
            preprocess(vcf_file=vcf_file, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file, anc_allele_file=None,
                       feature_config=kwargs["feature_config"], is_phased=kwargs["is_phased"],
                       ploidy=kwargs["ploidy"], output_dir=f'{kwargs["output_dir"]}/{rep}',
                       output_prefix=f'{kwargs["output_prefix"]}.{rep}', win_len=kwargs["seq_len"], win_step=kwargs["seq_len"], thread=1)

            feature_file = f'{kwargs["output_dir"]}/{rep}/{kwargs["output_prefix"]}.{rep}.features'
            feature_df = pd.read_csv(feature_file, sep="\t")
            feature_df['rep'] = rep
            feature_df.to_csv(feature_file, sep="\t", index=False)
            _label(ind_file=tgt_ind_file, truth_tract_file=bed_file, output=label_file,
                   is_phased=kwargs["is_phased"], ploidy=kwargs["ploidy"], rep=rep,
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

    Returns:
        tracts pandas.DataFrame: Ground truth introgressed fragments from a given source population to a given target populations.
    """
    tracts = pd.DataFrame(columns=['chrom', 'start', 'end', 'sample'])

    src_id = [p.id for p in ts.populations() if p.metadata['name']==src_id][0]
    tgt_id = [p.id for p in ts.populations() if p.metadata['name']==tgt_id][0]

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
                for n in ts.samples(tgt_id):
                    if t.is_descendant(n, m.node):
                        left = m.left if m.left > t.interval.left else t.interval.left
                        right = m.right if m.right < t.interval.right else t.interval.right
                        tracts.loc[len(tracts.index)] = [1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}']

    tracts = tracts.sort_values(by=['sample', 'chrom', 'start', 'end'])

    return tracts


def _label(ind_file, is_phased, ploidy, truth_tract_file, seq_len, intro_prop, not_intro_prop, rep, output):
    """
    """
    label_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'sample'])
    with open(ind_file, 'r') as f:
        if is_phased is True:
            for line in f:
                sample = line.rstrip()
                for p in range(ploidy):
                    label_df.loc[len(label_df.index)] = [1, 0, seq_len, f'{sample}_{p+1}']
        else:
            for line in f:
                sample = line.rstrip()
                label_df.loc[len(label_df.index)] = [1, 0, seq_len, sample]

    try:
        truth_tract_df = pd.read_csv(truth_tract_file, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        label_df['label'] = 0.0
    else:
        truth_tract_df.columns = ['chrom', 'start', 'end', 'sample']
        truth_tract_df['len'] = truth_tract_df['end'] - truth_tract_df['start']
        truth_tract_df = truth_tract_df.groupby(by=['sample'])['len'].sum().reset_index()
        truth_tract_df['prop'] = truth_tract_df['len'] / seq_len
        truth_tract_df['label'] = truth_tract_df.apply(lambda row: _add_label(row, intro_prop, not_intro_prop), axis=1)
        label_df = label_df.merge(truth_tract_df.drop(columns=['len', 'prop']),
                                  left_on=['sample'], right_on=['sample'], how='left').fillna(0)
    finally:
        label_df['rep'] = rep
        label_df['label'] = label_df['label'].astype('int8')
        label_df.to_csv(output, sep="\t", index=False)


def _add_label(row, intro_prop, not_intro_prop):
    """
    """
    if row['prop'] > intro_prop: return 1.0
    elif row['prop'] < not_intro_prop: return 0.0
    else: return -1.0


if __name__ == '__main__':
    simulate(demo_model_file="./examples/models/ArchIE_3D19.yaml", nrep=1, nref=50, ntgt=50, 
             ref_id='Ref', tgt_id='Tgt', src_id='Ghost', ploidy=2, seq_len=50000, mut_rate=1.25e-8, rec_rate=1e-8, thread=2,
             feature_config=None, is_phased=True, intro_prop=0.7, not_intro_prop=0.3, keep_sim_data=True,
             output_prefix='test', output_dir='./sstar/test7', seed=555)
