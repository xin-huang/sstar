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
from sstar.utils import multiprocessing_manager


def simulate(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, ploidy, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed):
    """
    """
    multiprocessing_manager(worker_func=_simulation_worker, nrep=nrep, thread=thread,
                            demo_model_file=demo_model_file, nref=nref, ntgt=ntgt, 
                            ref_id=ref_id, tgt_id=tgt_id, src_id=src_id, ploidy=ploidy,
                            seq_len=seq_len, mut_rate=mut_rate, rec_rate=rec_rate, 
                            output_prefix=output_prefix, output_dir=output_dir, seed=seed)


def _simulation_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        rep = in_queue.get()

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
            # Tree-sequences are sorted by the left ends of the intervals.
            # Can skip those tree-sequences are not overlapped with the interval of i.
            if i.left > t.interval.right: continue
            if i.right <= t.interval.left: break 
            for n in tracts.keys():
                if t.is_descendant(n, i.node): tracts[n].append([1, int(i.left), int(i.right), f'hap_{int(n%2)}', f'tsk_{ts.node(n).individual}'])

    return tracts


if __name__ == '__main__':
    simulate(demo_model_file="./examples/models/BonoboGhost_4K19.yaml", nrep=1000, nref=50, ntgt=50, 
             ref_id='Western', tgt_id='Bonobo', src_id='Ghost', ploidy=2, seq_len=50000, mut_rate=1e-8, rec_rate=1e-8, thread=2, 
             output_prefix='test', output_dir='./sstar/test', seed=913)
