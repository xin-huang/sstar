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
from multiprocessing import Process, Queue


def train(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed=None, train_archie=False):
    """
    """
    demo_graph = demes.load(demo_model_file)
    demography = msprime.Demography.from_demes(demo_graph)
    samples = [
        msprime.SampleSet(nref, ploidy=2, population=ref_id),
        msprime.SampleSet(ntgt, ploidy=2, population=tgt_id),
    ]

    # simulate data
    _manager(nrep, demography, samples, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed)

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


def _manager(nrep, demography, samples, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed):
    """
    """
    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    res = []
    in_queue, out_queue = Queue(), Queue()
    workers = [Process(target=_simulation, args=(in_queue, out_queue, demography, samples, seq_len, mut_rate, rec_rate, output_prefix, output_dir, seed)) for i in range(thread)]

    for i in range(nrep): in_queue.put(i)

    try:
        for w in workers: w.start()
        for i in range(nrep):
            item = out_queue.get()
            if item != '': res.append(item)
        for w in workers: w.terminate()
    finally:
        for w in workers: w.join()


def _simulation(in_queue, out_queue, demography, samples, seq_len, mut_rate, rec_rate, output_prefix, output_dir, seed):
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
        ts.dump(output_dir+'/'+output_prefix+f'{rep}.ts')
        out_queue.put(rep)


def _get_true_tracts():
    """
    """
    pass


if __name__ == '__main__':
    train("./examples/models/BonoboGhost_4K19_no_introgression.yaml", 10, 10, 10, 'Western', 'Bonobo', 10000, 1e-8, 1e-8, 2, 'test', './examples')
