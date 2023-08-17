import sstar
import utils
import os
import demes
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import nbinom

import train
import infer
import preprocess

import argparse, os, sys, signal
import shutil


def main(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed, folder_partitions, create_testdirs = False):
    
    #this variable determines how many predicition-directories are created within one folder
    nrep_per_folder = int(nrep / folder_partitions)

    train_df_list = []
    for i in range(folder_partitions):
        curr_output_dir = output_dir + str(i)

        preprocess.store_global_parameters(demo_model_file, nrep_per_folder, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, curr_output_dir)
        #create a training folder for the training set
        if not os.path.exists(curr_output_dir):
            os.makedirs(curr_output_dir)

        train._simulation_manager(demo_model_file, nrep_per_folder, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, curr_output_dir, seed)

        new_train_df = train._train_archie(demo_model_file, nrep_per_folder, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, curr_output_dir, drop_dynamic_cols=False)
        
        train_df_list.append(new_train_df)

        #after appending the features to the dataframe, the folder continaing training examples is deleted
        shutil.rmtree(curr_output_dir)

    #the full train dataframe
    train_df = pd.concat(train_df_list)

    #create reduced dfs for training on reduced data sets
    train_df_reduced = train_df.copy()
    train_df_no_kurtosis = train_df.copy()
    train_df_no_paired_dist = train_df.copy()
    train_df_target_full_reduced = train_df.copy()

    #drop_dynamic_cols indicate whether non-fixed size features should be dropped etc
    dynamic_cols = [col for col in train_df.columns if ('-ton' in col or col.startswith("pairwised_dist"))]
    no_kurtosis_cols = [col for col in train_df.columns if ('kurtosis_pairwised_dist' in col or col.startswith("pairwised_dist")) ]
    no_paired_cols = [col for col in train_df.columns if (col.startswith("pairwised_dist"))]
    full_reduced_cols = [col for col in train_df.columns if ('-ton' in col or 'pairwised_dist' in col )]

    train_df_reduced.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')
    train_df_no_kurtosis.drop(no_kurtosis_cols, axis=1, inplace = True, errors='ignore')
    train_df_no_paired_dist.drop(no_paired_cols, axis=1, inplace = True, errors='ignore')
    train_df_target_full_reduced.drop(full_reduced_cols, axis=1, inplace = True, errors='ignore')

    #reduced dataframes
    train_df.to_csv(str(demo_model_file) + "_nref" + str(nref) + "_ntgt" + str(ntgt) + "_finalfeaturefile.csv")
    train_df_reduced.to_csv(str(demo_model_file) + "_nref" + str(nref) + "_ntgt" + str(ntgt) + "_finalfeaturefile_fixed.csv")
    train_df_no_kurtosis.to_csv(str(demo_model_file) + "_nref" + str(nref) + "_ntgt" + str(ntgt) + "_finalfeaturefile_nokurtosis.csv")
    train_df_no_paired_dist.to_csv(str(demo_model_file) + "_nref" + str(nref) + "_ntgt" + str(ntgt) + "_finalfeaturefile_nopaired.csv")
    train_df_target_full_reduced.to_csv(str(demo_model_file) + "_nref" + str(nref) + "_ntgt" + str(ntgt) + "_finalfeaturefile_tgtfullreduced.csv")

    #names for models
    scikit_file = output_prefix + ".scikit.pickle"
    statsmodels_file = output_prefix + ".statsmodels.pickle"
    scikit_file_reduced = "fixed_" + output_prefix + ".scikit.pickle"
    statsmodels_file_reduced = "fixed_" + output_prefix + ".statsmodels.pickle"
    scikit_file_no_kurtosis = "nokurt_" + output_prefix + ".scikit.pickle"
    statsmodels_file_no_kurtosis = "nokurt_" + output_prefix + ".statsmodels.pickle"
    scikit_file_no_paired_dist = "nopaired_" + output_prefix + ".scikit.pickle"
    statsmodels_file_no_paired_dist = "nopaired_" + output_prefix + ".statsmodels.pickle"
    scikit_file_full_reduced = "fullreduced_" + output_prefix + ".scikit.pickle"
    statsmodels_file_full_reduced = "fullreduced_" + output_prefix + ".statsmodels.pickle"
    scikit_file_reduced = "fixed_" + output_prefix + ".scikit.pickle"
    statsmodels_file_reduced = "fixed_" + output_prefix + ".statsmodels.pickle"

    #call training functions

    train.train_statsmodels(train_df, statsmodels_file)
    train.train_scikit(train_df, scikit_file)

    train.train_statsmodels(train_df_reduced, statsmodels_file_reduced)
    train.train_scikit(train_df_reduced, scikit_file_reduced)

    train.train_statsmodels(train_df_no_kurtosis, statsmodels_file_no_kurtosis)
    train.train_scikit(train_df_no_kurtosis, scikit_file_no_kurtosis)

    train.train_statsmodels(train_df_no_paired_dist, statsmodels_file_no_paired_dist)
    train.train_scikit(train_df_no_paired_dist, scikit_file_no_paired_dist)

    train.train_statsmodels(train_df_target_full_reduced, statsmodels_file_full_reduced)
    train.train_scikit(train_df_target_full_reduced, scikit_file_full_reduced)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--demo_model_file',type=str, required=True)
    parser.add_argument('--nrep', type=int, required=True)
    parser.add_argument('--nref', type=int, required=True)
    parser.add_argument('--ntgt',type=int,  required=True)
    parser.add_argument('--ref_id',type=str, required=True)
    parser.add_argument('--tgt_id', type=str, required=True)
    parser.add_argument('--src_id', type=str, required=True)
    parser.add_argument('--seq_len', type=int, required=True)
    parser.add_argument('--mut_rate',type=float, required=True)
    parser.add_argument('--rec_rate',type=float, required=True)
    parser.add_argument('--thread',type=int, required=True)
    parser.add_argument('--output_prefix',type=str, required=True)
    parser.add_argument('--output_dir',type=str, required=True)
    parser.add_argument('--seed',required=True)

    parser.add_argument('--folder_partitions',type=int,required=True)
    
    args = parser.parse_args()
    demo_model_file = args.demo_model_file
    nrep = args.nrep
    nref = args.nref
    ntgt = args.ntgt
    ref_id = args.ref_id
    tgt_id = args.tgt_id
    src_id = args.src_id
    seq_len = args.seq_len
    mut_rate = args.mut_rate
    rec_rate = args.rec_rate
    thread = args.thread
    output_prefix = args.output_prefix
    output_dir = args.output_dir
    if args.seed == "None":
        seed = None
    else:
        seed = int(args.seed)

    folder_partitions = args.folder_partitions
    
    main(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed, folder_partitions)
