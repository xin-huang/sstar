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

#def main(args=None):
def main(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed, model_name):

    ref_ind_file = str(demo_model_file) + "_new_sim" + "_nref" + str(nref) + "_ntgt" + str(ntgt) + ".ref.ind.list"
    tgt_ind_file  = str(demo_model_file) + "_new_sim" + "_nref" + str(nref) + "_ntgt" + str(ntgt) + ".tgt.ind.list"

    scikitfile = output_prefix + ".scikit.pickle"
    statsmodelsfile = output_prefix + ".statsmodels.pickle"

    #get all folders for prediction
    final_folders = infer.get_all_folders(model_name, os.path.join("nref_" + str(nref), "ntgt_" + str(ntgt)))

    sample_name = "nref_" + str(nref) + "_ntgt_" + str(ntgt)

    #without ref_ and tgt_ind_file (are created within infer.predict_introgression_folders)
    #infer.predict_introgression_folders(nrep, nref, ntgt, seq_len, thread, output_prefix+ "test", final_folders, statsmodel=statsmodelsfile, scikitmodel=scikitfile, sample_name=sample_name, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file, model_name=model_name, drop_dynamic_cols=False, evaluate=False, simulated=True, average_for_inference=False, compute_cutoffs=True, win_step_50k = False)
    
    infer.predict_introgression_folders(nrep, nref, ntgt, seq_len, thread, output_prefix+ "test", final_folders, statsmodel=statsmodelsfile, scikitmodel=scikitfile, sample_name=sample_name, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file, model_name=model_name, drop_dynamic_cols=False, evaluate=False, simulated=True, average_for_inference=False, compute_cutoffs=True, win_step_50k = False)


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
    parser.add_argument('--model_name',type=str,required=True)

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
    model_name=args.model_name
    
    main(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed,model_name)
