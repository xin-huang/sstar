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
import numpy as np
import pandas as pd
import scipy
import yaml
from multiprocessing import Process, Queue
from sstar.stats import *
from sstar.utils import read_data, filter_data, create_windows, multiprocessing_manager


def preprocess(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, feature_config, is_phased, 
               ploidy, output_dir, output_prefix, win_len, win_step, thread):
    """
    Description:
        Processes genotype data.

    Arguments:
        vcf_file str: Name of the VCF file containing genotype data.
        ref_ind_file str: Name of the file containing sample information from the reference population.
        tgt_ind_file str: Name of the file containing sample information from the target population.
        anc_allele_file str: Name of the file containing ancestral allele information.
        feature_config str: Name of the YAML file specifying what features should be used. 
        is_phased bool: True, genomes are phased; False, genomes are unphased.
        ploidy int: Ploidy of genomes.
        output_dir str: Directory storing the output files.
        output_prefix str: Prefix of the output files.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        thread int: Number of threads.
    """
    if feature_config is None:
        features = None
        output_genotypes = True
    else:
        with open(feature_config, 'r') as f:
            features = yaml.safe_load(f)

        features = features['features']
        output_genotypes = False

    ref_data, ref_samples, tgt_data, tgt_samples = read_data(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, is_phased)

    for c in tgt_data.keys():
        windows = create_windows(tgt_data[c]['POS'], c, win_step, win_len)

    res = multiprocessing_manager(worker_func=_preprocess_worker, nrep=len(windows), thread=thread, windows=windows, ref_data=ref_data, tgt_data=tgt_data, 
                                  features=features, is_phased=is_phased, ploidy=ploidy, output_genotypes=output_genotypes)

    # x[0]: the chromosome name in number
    # x[1]: the start of the window
    # x[2]: the end of the window
    res.sort(key=lambda x: (x[0], x[1], x[2]))

    header = _create_header(ref_samples, tgt_samples, features, is_phased, ploidy, output_genotypes)
    _output(res, tgt_samples, header, features, is_phased, ploidy, output_dir, output_prefix, output_genotypes)


def _preprocess_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        chr_name, start, end, ref_gts, tgt_gts, pos = in_queue.get()
        variants_not_in_ref = np.sum(ref_gts, axis=1) == 0
        sub_ref_gts = ref_gts[variants_not_in_ref]
        sub_tgt_gts = tgt_gts[variants_not_in_ref]
        sub_pos = pos[variants_not_in_ref]

        items = dict()

        if kwargs['output_genotypes'] is True:
            items['ref_gts'] = ref_gts
            items['tgt_gts'] = tgt_gts
        else:
            if 'number of total mutations' in kwargs['features'].keys():
                items['ttl_mut_nums'] = cal_mut_num(ref_gts, tgt_gts, mut_type='total')
            if 'number of private mutations' in kwargs['features'].keys():
                items['pvt_mut_nums'] = cal_mut_num(sub_ref_gts, sub_tgt_gts, mut_type='private')
            if 'individual allele frequency spectra' in kwargs['features'].keys():
                if kwargs['is_phased'] is True: ploidy = 1
                else: ploidy = kwargs['ploidy']
                spectra = cal_n_ton(tgt_gts, ploidy=ploidy)
                items['spectra'] = spectra
            if 'reference distances' in kwargs['features'].keys():
                ref_dists = cal_dist(ref_gts, tgt_gts)
                if 'all' in kwargs['features']['reference distances'].keys(): items['ref_dists'] = ref_dists
                if 'minimum' in kwargs['features']['reference distances'].keys(): items['min_ref_dists'] = np.min(ref_dists, axis=1)
                if 'maximum' in kwargs['features']['reference distances'].keys(): items['max_ref_dists'] = np.max(ref_dists, axis=1)
                if 'mean' in kwargs['features']['reference distances'].keys(): items['mean_ref_dists'] = np.mean(ref_dists, axis=1)
                if 'median' in kwargs['features']['reference distances'].keys(): items['median_ref_dists'] = np.median(ref_dists, axis=1)
                if 'variance' in kwargs['features']['reference distances'].keys(): items['var_ref_dists'] = np.var(ref_dists, axis=1)
                if 'skew' in kwargs['features']['reference distances'].keys(): items['skew_ref_dists'] = scipy.stats.skew(ref_dists, axis=1)
                if 'kurtosis' in kwargs['features']['reference distances'].keys(): items['kurtosis_ref_dists'] = scipy.stats.kurtosis(ref_dists, axis=1)
            if 'target distances' in kwargs['features'].keys():
                tgt_dists = cal_dist(tgt_gts, tgt_gts)
                if 'all' in kwargs['features']['target distances'].keys(): items['tgt_dists'] = tgt_dists
                if 'minimum' in kwargs['features']['target distances'].keys(): items['min_tgt_dists'] = np.min(tgt_dists, axis=1)
                if 'maximum' in kwargs['features']['target distances'].keys(): items['max_tgt_dists'] = np.max(tgt_dists, axis=1)
                if 'mean' in kwargs['features']['target distances'].keys(): items['mean_tgt_dists'] = np.mean(tgt_dists, axis=1)
                if 'median' in kwargs['features']['target distances'].keys(): items['median_tgt_dists'] = np.median(tgt_dists, axis=1)
                if 'variance' in kwargs['features']['target distances'].keys(): items['var_tgt_dists'] = np.var(tgt_dists, axis=1)
                if 'skew' in kwargs['features']['target distances'].keys(): items['skew_tgt_dists'] = scipy.stats.skew(tgt_dists, axis=1)
                if 'kurtosis' in kwargs['features']['target distances'].keys(): items['kurtosis_tgt_dists'] = scipy.stats.kurtosis(tgt_dists, axis=1)
            if 'sstar' in kwargs['features'].keys():
                sstar_scores, sstar_snp_nums, haplotypes = cal_sstar(sub_tgt_gts, sub_pos, 
                                                                     method=kwargs['features']['sstar']['genotype distance'], 
                                                                     match_bonus=kwargs['features']['sstar']['match bonus'],
                                                                     max_mismatch=kwargs['features']['sstar']['max mismatch'],
                                                                     mismatch_penalty=kwargs['features']['sstar']['mismatch penalty'])
                items['sstar'] = sstar_scores

        out_queue.put((chr_name, start, end, items))


def _create_header(ref_samples, tgt_samples, features, is_phased, ploidy, output_genotypes):
    """
    """
    if output_genotypes is True:
        if is_phased is True:
            haps = []
            for s in ref_samples:
                for i in range(ploidy):
                    haps.append(f'{s}_{i+1}')
            for s in tgt_samples:
                for i in range(ploidy):
                    haps.append(f'{s}_{i+1}')
            header = "\t".join(haps)
        else: header = "\t".join(ref_samples) + "\t" + "\t".join(tgt_samples)
    else:
        header = "chrom\tstart\tend\tsample"
        if 'sstar' in features.keys(): header += "\tS*_score"
        if 'number of total mutations' in features.keys(): header += "\ttotal_SNP_num"
        if 'number of private mutations' in features.keys(): header += "\tprivate_SNP_num"
        if 'individual allele frequency spectra' in features.keys():
            nsample = len(tgt_samples)*ploidy
            for i in range(nsample+1): header += f'\t{i}_ton'
        if 'reference distances' in features.keys():
            if 'minimum' in features['reference distances'].keys(): header += "\tmin_ref_dist"
            if 'maximum' in features['reference distances'].keys(): header += "\tmax_ref_dist"
            if 'mean' in features['reference distances'].keys(): header += "\tmean_ref_dist"
            if 'median' in features['reference distances'].keys(): header += "\tmedian_ref_dist"
            if 'variance' in features['reference distances'].keys(): header += "\tvar_ref_dist"
            if 'skew' in features['reference distances'].keys(): header += "\tskew_ref_dist"
            if 'kurtosis' in features['reference distances'].keys(): header += "\tkurtosis_ref_dist"
            if 'all' in features['reference distances'].keys(): 
                if is_phased is True: 
                    for s in ref_samples:
                        for i in range(ploidy):
                            header += f'\tref_dist_{s}_{i+1}'
                else:
                    for s in ref_samples:
                        header += f'\tref_dist_{s}'
        if 'target distances' in features.keys():
            if 'minimum' in features['target distances'].keys(): header += "\tmin_tgt_dist"
            if 'maximum' in features['target distances'].keys(): header += "\tmax_tgt_dist"
            if 'mean' in features['target distances'].keys(): header += "\tmean_tgt_dist"
            if 'median' in features['target distances'].keys(): header += "\tmedian_tgt_dist"
            if 'variance' in features['target distances'].keys(): header += "\tvar_tgt_dist"
            if 'skew' in features['target distances'].keys(): header += "\tskew_tgt_dist"
            if 'kurtosis' in features['target distances'].keys(): header += "\tkurtosis_tgt_dist"
            if 'all' in features['target distances'].keys(): 
                if is_phased is True: 
                    for s in tgt_samples:
                        for i in range(ploidy):
                            header += f'\ttgt_dist_{s}_{i+1}'
                else: 
                    for s in tgt_samples:
                        header += f'\ttgt_dist_{s}'

    return header


def _output(res, tgt_samples, header, features, is_phased, ploidy, output_dir, output_prefix, output_genotypes):
    """
    """
    os.makedirs(output_dir, exist_ok=True)
    if output_genotypes is True:
        os.makedirs(f'{output_dir}/genotypes', exist_ok=True)
        for r in res:
            chrom = r[0]
            start = r[1]
            end = r[2]
            items = r[3]
            os.makedirs(f'{output_dir}/genotypes/{chrom}', exist_ok=True)
            if is_phased is True: output_file = f'{output_dir}/genotypes/{chrom}/{output_prefix}.{chrom}.{start}-{end}.phased.genotypes'
            else: output_file = f'{output_dir}/genotypes/{chrom}/{output_prefix}.{chrom}.{start}-{end}.unphased.genotypes'
            with open(output_file, 'w') as f:
                f.write(f'{header}\n')
                for i in range(len(items['ref_gts'])):
                    ref_gts = "\t".join(items['ref_gts'][i].astype(str))
                    tgt_gts = "\t".join(items['tgt_gts'][i].astype(str))
                    f.write(f'{ref_gts}\t{tgt_gts}\n')
    else:
        output_file = f'{output_dir}/{output_prefix}.features'
        if is_phased is not True: ploidy = 1
        with open(output_file, 'w') as f:
            f.write(f'{header}\n')
            for r in res:
                chrom = r[0]
                start = r[1]
                end = r[2]
                items = r[3]
                for i in range(len(tgt_samples)*ploidy):
                    if ploidy != 1: sample = f'{tgt_samples[int(i/ploidy)]}_{i%ploidy+1}'
                    else: sample = tgt_samples[i]
                    out = ''
                    if 'sstar' in features.keys(): out += f'{items["sstar"][i]}'
                    if 'number of total mutations' in features.keys(): out += f'\t{items["ttl_mut_nums"][i]}'
                    if 'number of private mutations' in features.keys(): out += f'\t{items["pvt_mut_nums"][i]}'
                    if 'individual allele frequency spectra' in features.keys():
                        spectra = "\t".join(items["spectra"][i].astype(str))
                        out += f'\t{spectra}'
                    if 'reference distances' in features.keys():
                        if 'minimum' in features['reference distances'].keys(): out += f'\t{items["min_ref_dists"][i]}'
                        if 'maximum' in features['reference distances'].keys(): out += f'\t{items["max_ref_dists"][i]}'
                        if 'mean' in features['reference distances'].keys(): out += f'\t{items["mean_ref_dists"][i]}'
                        if 'median' in features['reference distances'].keys(): out += f'\t{items["median_ref_dists"][i]}'
                        if 'variance' in features['reference distances'].keys(): out += f'\t{items["var_ref_dists"][i]}'
                        if 'skew' in features['reference distances'].keys(): out += f'\t{items["skew_ref_dists"][i]}'
                        if 'kurtosis' in features['reference distances'].keys(): out += f'\t{items["kurtosis_ref_dists"][i]}'
                        if 'all' in features['reference distances'].keys(): 
                            dists = "\t".join(items["ref_dists"][i].astype(str))
                            out += f'\t{dists}'
                    if 'target distances' in features.keys():
                        if 'minimum' in features['target distances'].keys(): out += f'\t{items["min_tgt_dists"][i]}'
                        if 'maximum' in features['target distances'].keys(): out += f'\t{items["max_tgt_dists"][i]}'
                        if 'mean' in features['target distances'].keys(): out += f'\t{items["mean_tgt_dists"][i]}'
                        if 'median' in features['target distances'].keys(): out += f'\t{items["median_tgt_dists"][i]}'
                        if 'variance' in features['target distances'].keys(): out += f'\t{items["var_tgt_dists"][i]}'
                        if 'skew' in features['target distances'].keys(): out += f'\t{items["skew_tgt_dists"][i]}'
                        if 'kurtosis' in features['target distances'].keys(): out += f'\t{items["kurtosis_tgt_dists"][i]}'
                        if 'all' in features['target distances'].keys(): 
                            dists = "\t".join(items["tgt_dists"][i].astype(str))
                            out += f'\t{dists}'
                    f.write(f'{chrom}\t{start}\t{end}\t{sample}\t{out}\n')


def training_preprocess(input_dir, input_prefix, nrep, feature_config, is_phased,
                        ploidy, output_dir, output_prefix, seq_len, thread, intro_prop, not_intro_prop):
    """
    """
    res = multiprocessing_manager(worker_func=_training_preprocess_worker, nrep=nrep, thread=thread,
                                  feature_config=feature_config, is_phased=is_phased, ploidy=ploidy,
                                  seq_len=seq_len, intro_prop=intro_prop, not_intro_prop=not_intro_prop,
                                  input_dir=input_dir, input_prefix=input_prefix, output_dir=output_dir, output_prefix=output_prefix)

    if feature_config is not None:
        feature_df = pd.DataFrame()
        for i in range(nrep):
            df = pd.read_csv(f'{output_dir}/{i}/{output_prefix}.{i}.labeled.features', sep="\t")
            feature_df = pd.concat([feature_df, df])

        feature_df.to_csv(f'{output_dir}/{output_prefix}.all.labeled.features', sep="\t", index=False)


def _training_preprocess_worker(in_queue, out_queue, **kwargs):
    """
    """
    while True:
        rep, seed = in_queue.get()
        vcf = f'{kwargs["input_dir"]}/{rep}/{kwargs["input_prefix"]}.{rep}.vcf'
        ref = f'{kwargs["input_dir"]}/{rep}/{kwargs["input_prefix"]}.{rep}.ref.ind.list'
        tgt = f'{kwargs["input_dir"]}/{rep}/{kwargs["input_prefix"]}.{rep}.tgt.ind.list'
        output_dir = f'{kwargs["output_dir"]}/{rep}'
        output_prefix = f'{kwargs["output_prefix"]}.{rep}'

        preprocess(vcf_file=vcf, ref_ind_file=ref, tgt_ind_file=tgt, anc_allele_file=None,
                   feature_config=kwargs["feature_config"], is_phased=kwargs["is_phased"], 
                   ploidy=kwargs["ploidy"], output_dir=output_dir,
                   output_prefix=output_prefix, win_len=kwargs["seq_len"], win_step=kwargs["seq_len"], thread=1)

        feature_file = f'{output_dir}/{output_prefix}.features'
        truth_tract_file = f'{output_dir}/{output_prefix}.truth.tracts.bed'
        output = f'{output_dir}/{output_prefix}.labeled.features'

        _label(feature_file=feature_file, truth_tract_file=truth_tract_file, output=output,
               seq_len=kwargs["seq_len"], intro_prop=kwargs["intro_prop"], not_intro_prop=kwargs["not_intro_prop"])

        out_queue.put(rep)


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
    #preprocess(vcf_file="examples/data/real_data/sstar.example.biallelic.snps.vcf.gz", ref_ind_file="examples/data/ind_list/ref.ind.list", tgt_ind_file="examples/data/ind_list/tgt.ind.list", anc_allele_file=None, feature_config="examples/features/sstar.features.yaml", is_phased=False, ploidy=2, output_dir="./test_data", output_prefix="test.sstar", win_len=50000, win_step=10000, thread=1)
    training_preprocess(input_dir="./sstar/test", input_prefix="test", nrep=1000, feature_config="./examples/features/archie.features.yaml", is_phased=True,
                        ploidy=2, output_dir="./sstar/test", output_prefix="test", seq_len=50000, thread=2, archaic_prop=0.7, not_archaic_prop=0.3)
