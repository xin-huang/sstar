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

import os
import demes
import glob
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Process, Queue
from scipy.stats import norm
from scipy.stats import nbinom


def get_quantile(model, ms_dir, N0, nsamp, nreps, ref_index, ref_size, tgt_index, tgt_size, mut_rate, rec_rate, seq_len, snp_num_range, output_dir, thread, all_ind_geno_dist, seeds):
    """
    """
    if seeds is not None: np.random.seed(np.sum(seeds))
    output_dir = os.path.abspath(output_dir)
    if os.path.exists(output_dir) is False: subprocess.call(['mkdir', output_dir])
    _generate_mut_rec_combination(N0, nreps, mut_rate, rec_rate, seq_len, output_dir)
    _run_ms_simulation(model, ms_dir, N0, nsamp, nreps, ref_index, ref_size, tgt_index, tgt_size, seq_len, snp_num_range, output_dir, thread, all_ind_geno_dist, seeds)
    _summary(output_dir, rec_rate)


def _generate_mut_rec_combination(N0, nreps, mut_rate, rec_rate, seq_len, output_dir):
    """
    """
    scaled_mut_rate = 4*N0*mut_rate*seq_len
    scaled_rec_rate = 4*N0*rec_rate*seq_len
    mut_rate_list = norm.rvs(loc=scaled_mut_rate, scale=0.233, size=nreps)
    rec_rate_list = nbinom.rvs(n=0.5, p=0.5/(0.5+scaled_rec_rate), size=nreps)

    rates = f'{output_dir}/rates.combination'
    with open(rates, 'w') as o:
        for i in range(len(mut_rate_list)):
            if mut_rate_list[i] < 0.001: mut_rate_list[i] = 0.001
            if rec_rate_list[i] < 0.001: rec_rate_list[i] = 0.001
            mut_rate = mut_rate_list[i]
            rec_rate = rec_rate_list[i]
            o.write(f'{mut_rate}\t{rec_rate}\n')


def _run_ms_simulation(model, ms_dir, N0, nsamp, nreps, ref_index, ref_size, tgt_index, tgt_size, seq_len, snp_num_range, output_dir, thread, all_ind_geno_dist, seeds):
    """
    """
    graph = demes.load(model)
    samples = np.zeros(len(graph.demes))
    samples[ref_index-1] = ref_size
    samples[tgt_index-1] = tgt_size
    ms_params = demes.to_ms(graph, N0=N0, samples=samples)
    snp_num_list = np.arange(snp_num_range[0], snp_num_range[1]+snp_num_range[2], snp_num_range[2])
    
    ms_exec = os.path.abspath(ms_dir) + '/ms'
    rates = f'{output_dir}/rates.combination'
    ref_list = f'{output_dir}/sim.ref.list'
    tgt_list = f'{output_dir}/sim.tgt.list'
    ref_size = int(ref_size / 2)
    tgt_size = int(tgt_size / 2)

    if ref_index == tgt_index: raise Exception('The reference population should be different from the target population.')
    elif ref_index < tgt_index:
        with open(ref_list, 'w') as o:
            for i in range(ref_size):
                o.write(f'ms_{i}\n')
        with open(tgt_list, 'w') as o:
            for i in range(tgt_size):
                o.write(f'ms_{i+ref_size}\n')
    else:
        with open(ref_list, 'w') as o:
            for i in range(ref_size):
                o.write(f'ms_{i+tgt_size}\n')
        with open(tgt_list, 'w') as o:
            for i in range(tgt_size):
                o.write(f'ms_{i}\n')
       
    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    in_queue, out_queue = Queue(), Queue()
    workers = [Process(target=_run_ms_simulation_worker, args=(in_queue, out_queue, output_dir, rates, ms_exec, nsamp, nreps, seq_len, ms_params, ref_list, tgt_list, all_ind_geno_dist, seeds)) for ii in range(thread)]
 
    for snp_num in snp_num_list:
        in_queue.put(snp_num)

    try:
        for worker in workers:
            worker.start()
        for snp_num in snp_num_list:
            out_queue.get()
        for worker in workers:
            worker.terminate()
    finally:
        for worker in workers:
            worker.join()


def _run_ms_simulation_worker(in_queue, out_queue, output_dir, rates, ms_exec, nsamp, nreps, seq_len, ms_params, ref_list, tgt_list, all_ind_geno_dist, seeds):
    """
    """
    while True:
        snp_num = in_queue.get()
        output_subdir = f'{output_dir}/{snp_num}'
        output_ms = f'{output_subdir}/sim.ms'
        output_vcf = f'{output_subdir}/sim.vcf'
        output_score = f'{output_subdir}/sim.score'
        output_quantile = f'{output_subdir}/sim.quantile'
        ms_script = f'{output_subdir}/run_ms.sh'
        cmd = " ".join(['cat', rates, '|', ms_exec, str(nsamp), str(nreps), '-seeds', " ".join([str(s) for s in seeds]), '-t', 'tbs', '-r', 'tbs', str(seq_len), '-s', str(snp_num), ms_params, '>', output_ms])
        
        if os.path.exists(output_subdir) is False: subprocess.call(['mkdir', output_subdir])
        with open(ms_script, 'w') as o:
            o.write(cmd+"\n")
        subprocess.call(['bash', ms_script])
        _ms2vcf(output_ms, output_vcf, nsamp, seq_len)
        if all_ind_geno_dist is True: subprocess.call(['sstar', 'score', '--vcf', output_vcf, '--ref', ref_list, '--tgt', tgt_list, '--output', output_score, '--win-len', str(seq_len), '--win-step', str(seq_len), '--thread', '1', '--all-ind-geno-dist'])
        else: subprocess.call(['sstar', 'score', '--vcf', output_vcf, '--ref', ref_list, '--tgt', tgt_list, '--output', output_score, '--win-len', str(seq_len), '--win-step', str(seq_len), '--thread', '1'])
        _cal_quantile(output_score, output_quantile, snp_num)
        out_queue.put('Finished')


def _ms2vcf(ms_file, vcf_file, nsamp, seq_len, ploidy=2):
    """
    Description:
        Converts ms output files into the VCF format.

    Arguments:
        ms_file str: Name of the ms file (input).
        vcf_file str: Name of the VCF file (output).
        nsamp int: Number of haploid genomes.
        seq_len int: Sequence length.
        ploidy int: Ploidy of each individual.
    """
    data = []
    i = -1
    header = "##fileformat=VCFv4.2\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(['ms_' + str(i) for i in range(int(nsamp/ploidy))])

    with open(ms_file, 'r') as f:
        f.readline()
        f.readline()
        for l in f.readlines():
            if l.startswith('//'):
                i += 1
                data.append({})
                data[i]['pos'] = []
                data[i]['geno'] = []
            elif l.startswith('positions'):
                data[i]['pos'] = l.rstrip().split(" ")[1:]
            elif l.startswith('0') or l.startswith('1'):
                data[i]['geno'].append(l.rstrip())

    shift = 0
    with open(vcf_file, 'w') as o:
        o.write(header+"\n")
        for i in range(len(data)):
            for j in range(len(data[i]['pos'])):
                pos = int(seq_len * float(data[i]['pos'][j])) + shift
                genotypes = "".join([data[i]['geno'][k][j] for k in range(len(data[i]['geno']))])
                genotypes = "\t".join([a+'|'+b for a,b in zip(genotypes[0::ploidy],genotypes[1::ploidy])])
                o.write(f"1\t{pos}\t.\tA\tT\t100\tPASS\t.\tGT\t{genotypes}\n")
            shift += seq_len


def _cal_quantile(in_file, out_file, snp_num):
    """
    """
    df = pd.read_csv(in_file, sep="\t").dropna()
    quantiles = np.arange(0.5,1,0.005)
    mean_df = df.groupby(['chrom', 'start', 'end'], as_index=False)['S*_score'].mean().dropna()
    scores = np.quantile(mean_df['S*_score'], quantiles)
    with open(out_file, 'w') as o:
        o.write('S*_score\tSNP_num\tquantile\n')
        for i in range(len(scores)):
            o.write(f'{scores[i]}\t{snp_num}\t{quantiles[i]}\n')


def _summary(output_dir, rec_rate):
    """
    """
    all_files = glob.glob(f'{output_dir}/*/*.quantile')
    li = []
    for filename in all_files:
        df = pd.read_csv(filename, sep='\t')
        li.append(df)

    df = pd.concat(li, ignore_index=True).round(3)
    df['rec_rate'] = rec_rate
    df.sort_values(by=['SNP_num', 'quantile']).to_csv(f'{output_dir}/quantile.summary.txt', sep='\t', index=False)
