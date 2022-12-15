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


import demes, msprime
import pandas as pd
from multiprocessing import Process, Queue


def train(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed=None, train_archie=False):
    """
    """
    # simulate data
    _simulation_manager(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed)

    if train_archie:
        _train_archie()
    else:
        _train_sstar()


def _train_archie():
    """
    """
    pass


def _train_sstar():
    pass


def _simulation_manager(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed):
    """
    """
    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    demo_graph = demes.load(demo_model_file)
    demography = msprime.Demography.from_demes(demo_graph)
    samples = [
        msprime.SampleSet(nref, ploidy=2, population=ref_id),
        msprime.SampleSet(ntgt, ploidy=2, population=tgt_id),
    ]

    res = []
    in_queue, out_queue = Queue(), Queue()
    workers = [Process(target=_simulation_worker, args=(in_queue, out_queue, demography, samples, tgt_id, src_id, seq_len, mut_rate, rec_rate, output_prefix, output_dir, seed)) for i in range(thread)]

    for i in range(nrep): in_queue.put(i)

    try:
        for w in workers: w.start()
        for i in range(nrep):
            item = out_queue.get()
            if item != '': res.append(item)
        for w in workers: w.terminate()
    finally:
        for w in workers: w.join()


def _simulation_worker(in_queue, out_queue, demography, samples, tgt_id, src_id, seq_len, mut_rate, rec_rate, output_prefix, output_dir, seed):
    """
    """
    while True:
        rep = in_queue.get()
        ts = msprime.sim_ancestry(
            recombination_rate=rec_rate,
            sequence_length=seq_len,
            samples=samples,
            demography=demography,
            record_migrations=True,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed)
        true_tracts = _get_true_tracts(ts, tgt_id, src_id)

        ts.dump(output_dir+'/'+output_prefix+f'{rep}.ts')
        with open(output_dir+'/'+output_prefix+f'{rep}.vcf', 'w') as o:
            ts.write_vcf(o)
       
        df = pd.DataFrame()
        for n in true_tracts.keys():
            true_tracts[n].sort(key=lambda x:(x[0], x[1], x[2]))
            df2 = pd.DataFrame(true_tracts[n], columns=['chr', 'start', 'end', 'hap', 'ind'])
            df = pd.concat([df, df2])

        df.drop_duplicates(keep='first').to_csv(output_dir+'/'+output_prefix+f'{rep}.true.tracts.bed', sep="\t", header=False, index=False)

        out_queue.put(rep)


def _get_true_tracts(ts, tgt_id, src_id):
    """
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


if __name__ == '__main__':
    train("./examples/models/BonoboGhost_4K19.yaml", 1000, 50, 50, 'Western', 'Bonobo', 'Ghost', 50000, 1e-8, 1e-8, 2, 'test', './sstar/test')
