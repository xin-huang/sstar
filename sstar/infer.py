
#from sstar import stats
#from sstar import preprocess
#from sstar import train
import stats
import preprocess
import train

import os
import allel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_all_folders(output_dir, ref_tgt_folder):
    res_dir = os.path.join("results", "simulated_data", output_dir, ref_tgt_folder)
    rep_folders = []
    for replicate, folder in enumerate(os.listdir(res_dir)):
        rep_folders.append(os.path.join(res_dir, folder))
    return rep_folders


def simulate_predict_introgression(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, statsmodel=None, scikitmodel=None):
    output_dir = output_dir + "test"
    output_prefix = output_prefix + "test"
    create_testdata(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir)

    predict_introgression(seq_len, thread,output_prefix, output_dir, statsmodel, scikitmodel)


def predict_introgression(seq_len, thread, output_prefix, output_dir, statsmodel=None, scikitmodel=None, drop_dynamic_cols=True, evaluate=True, simulated=True, average_for_inference=False, compute_cutoffs=True, ref_ind_file=None, tgt_ind_file=None, compute_statsmodel =False, plot_curves = False):    
    #set filenames for individuals, reference and target (needed for prprocess.process_data)

    if ref_ind_file == None:
        ref_ind_file = "new_sim.ref.ind.list"
    if tgt_ind_file == None:
        tgt_ind_file  = "new_sim.tgt.ind.list"


    anc_allele_file = None

    #set window length and stepsize
    win_len = 50000
    win_step = 10000
 
    #I think these parameters are NOT necessary for ArchIE - just retained for the signature of preprocess.process_data
    match_bonus = 1
    max_mismatch = 1
    mismatch_penalty = 1

    process_archie = True

    #tracts with a proportion between not_archaic and archaic are labeled as ambiguous (in _label)
    archaic_prop = 0.7
    not_archaic_prop = 0.3
    

    true_tracts = []
    true_tracts_infer = []

    features = []
    file_names = []
    replicate_counter = 0

    for replicate, file in enumerate(os.listdir(output_dir)):
        if file.endswith(".vcf") or file.endswith("biallelic.vcf.gz"):
        #if file.endswith(".vcf") or (file.endswith("vcf.gz") and "biallelic" not in file):  

            filename = os.path.splitext(file)[0]
            feature_file = os.path.splitext(file)[0]+'.features'

            #computing of statistics etc.
            preprocess.process_data(os.path.join(output_dir,file), ref_ind_file, tgt_ind_file, anc_allele_file, os.path.join(output_dir,feature_file), win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, process_archie)

            #load true tracts
            if file.endswith(".vcf"):
                true_tract = os.path.splitext(file)[0]+'.true.tracts.bed'
            else:
                true_tract = file.split(".")[0]+'.introgressed.tracts.bed'
            
            true_tract_data = pd.read_csv(os.path.join(output_dir, true_tract), sep="\t", header=None, names=['chr', 'start', 'end', 'hap', 'ind'])

            add_ind = None
            if true_tract_data is not None:
                true_tract_data["rep"] = replicate_counter
                if file.endswith(".vcf"):
                    true_tract_data["hap"] = true_tract_data["hap"].str.replace("hap_", "")
                else:
                    
                    add_ind = "tsk_10"
                    true_tract_data["hap"] = 1
                    true_tract_data["ind"] = add_ind


            true_tracts_infer.append ( true_tract_data )

            #label the true tracts according to the function from train
            true_tract_labeled = train._label(os.path.join(output_dir, true_tract), archaic_prop, not_archaic_prop, seq_len, add_ind)
            #true_tract_labeled = train._label(os.path.join(output_dir, true_tract), archaic_prop, not_archaic_prop, seq_len)
            if true_tract_labeled is not None:
                true_tract_labeled["rep"] = replicate_counter

            true_tracts.append(true_tract_labeled)

            #load the feature files with statistics created before
            feature = pd.read_csv(os.path.join(output_dir, feature_file), sep="\t")
            feature["rep"] = replicate_counter

            features.append(feature)
            file_names.append(filename)
            replicate_counter = replicate_counter + 1


    feature_df_labeleds = []

    for i, feature_df in enumerate(features):

        if compute_cutoffs == True:
            #label function from train is used
            feature_df_labeled = train.label_feature_df(feature_df, true_tracts_infer[i], discard_ambiguous=False, replicates=True)
            #feature_df_labeled = train.label_feature_df_archie(feature_df, true_tracts_infer[i], discard_ambiguous=False, replicates=True)

        else:

            feature_df_labeled = label_feature_df_archie_infer_only_label(feature_df, true_tracts[i], discard_ambiguous=False)


        if not os.path.exists(os.path.join(output_dir, "feature_dfs")):
            os.makedirs(os.path.join(output_dir, "feature_dfs"))

        feature_df_labeled.to_csv(os.path.join(output_dir, "feature_dfs", "features" + str(file_names[i]) + ".csv") )

        feature_df_labeleds.append(feature_df_labeled)

    #create one big test dataframe 

    test_df = pd.concat(feature_df_labeleds)

    #drop_dynamic_cols indicates whether all non-fixed size features should be dropped
    if drop_dynamic_cols == True:
        dynamic_cols = [col for col in test_df.columns if ('-ton' in col or col.startswith("pairwised_dist"))]

        test_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')

    #in the case of simulated data, we know the labels, but for inference we remove them for the test dataframe
    if simulated == True and average_for_inference == False:
        y_true = test_df["label"]
        test_df.drop(["label"], axis=1, inplace=True, errors='ignore')
    

    
    #statsmodel inference
    if statsmodel != None and compute_statsmodel==True:
        test_df_pred = test_df.copy()
        model_statsmodel = load_statsmodel(statsmodel)

        #only compute probabilities
        if compute_cutoffs == False:
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:       
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

            test_df_pred["probabilities"] = y_prob
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            new_dfs = split_dfs_sample_replicate(df_final)
            
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            
            if plot_curves == True:
                plot_cutoffs(recalls, precisions)


        #statsmodel inference
    if statsmodel != None and compute_statsmodel==True:
        test_df_pred = test_df.copy()
        model_statsmodel = load_statsmodel(statsmodel)

        #only compute probabilities
        if compute_cutoffs == False:
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:       
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

            test_df_pred["probabilities"] = y_prob
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            new_dfs = split_dfs_sample_replicate(df_final)
            
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            
            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / Statsmodel unphased")
    

    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit)
        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)

            if plot_curves == True:
                plot_cutoffs(recalls, precisions)
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            if evaluate == True:
                test_df_pred = test_df.copy()
                evaluate_scikit(test_df_pred, y_true, model_scikit)
       


    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit)
        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)

            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit unphased")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            #if evaluate == True:
            #    evaluate_scikit(test_df_pred, y_true, model_scikit)

def predict_introgression_folders_wo_parallel(seq_len, thread, output_prefix, output_dirs, statsmodel=None, scikitmodel=None, drop_dynamic_cols=True, evaluate=False, simulated=True, average_for_inference=False, compute_cutoffs=True, ref_ind_file=None, tgt_ind_file=None, model_name="archie", sample_name="sample1", compute_statsmodel = False, plot_curves=False):    
    #set filenames for individuals, reference and target (needed for prprocess.process_data)

    if ref_ind_file == None:
        ref_ind_file = "new_sim.ref.ind.list"
    if tgt_ind_file == None:
        tgt_ind_file  = "new_sim.tgt.ind.list"


    anc_allele_file = None

    #set window length and stepsize
    win_len = 50000
    win_step = 10000
 
    #I think these parameters are NOT necessary for ArchIE - just retained for the signature of preprocess.process_data
    match_bonus = 1
    max_mismatch = 1
    mismatch_penalty = 1

    process_archie = True

    #tracts with a proportion between not_archaic and archaic are labeled as ambiguous (in _label)
    archaic_prop = 0.7
    not_archaic_prop = 0.3
    

    true_tracts = []
    true_tracts_infer = []

    features = []
    file_names = []
    replicate_counter = 0

    feature_df_labeleds = []

    for output_dir in output_dirs:
        if os.path.isdir(output_dir):
            for replicate, file in enumerate(os.listdir(output_dir)):
            
                if file.endswith(".vcf") or file.endswith("biallelic.vcf.gz"):
                #if file.endswith(".vcf") or (file.endswith("vcf.gz") and "biallelic" not in file):  

                    filename = os.path.splitext(file)[0]
                    feature_file = os.path.splitext(file)[0]+'.features'

                    #computing of statistics etc.
                    preprocess.process_data(os.path.join(output_dir,file), ref_ind_file, tgt_ind_file, anc_allele_file, os.path.join(output_dir,feature_file), win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, process_archie)

                    #load true tracts
                    if file.endswith(".vcf"):
                        true_tract = os.path.splitext(file)[0]+'.true.tracts.bed'
                    else:
                        true_tract = file.split(".")[0]+'.introgressed.tracts.bed'
                    
                    true_tract_data = pd.read_csv(os.path.join(output_dir, true_tract), sep="\t", header=None, names=['chr', 'start', 'end', 'hap', 'ind'])


                    add_ind = None
                    if true_tract_data is not None:
                        true_tract_data["rep"] = replicate_counter
                        if file.endswith(".vcf"):
                            true_tract_data["hap"] = true_tract_data["hap"].str.replace("hap_", "")
                        else:
                        
                            add_ind = "tsk_10"
                            true_tract_data["hap"] = 1
                            true_tract_data["ind"] = add_ind


                    true_tracts_infer.append ( true_tract_data )

                    #label the true tracts according to the function from train
                    true_tract_labeled = train._label(os.path.join(output_dir, true_tract), archaic_prop, not_archaic_prop, seq_len, add_ind)

                    #label the true tracts according to the function from train
                    #true_tract_labeled = train._label(os.path.join(output_dir, true_tract), archaic_prop, not_archaic_prop, seq_len)
                    if true_tract_labeled is not None:
                        true_tract_labeled["rep"] = replicate_counter

                    true_tracts.append(true_tract_labeled)

                    #load the feature files with statistics created before
                    feature = pd.read_csv(os.path.join(output_dir, feature_file), sep="\t")
                    feature["rep"] = replicate_counter

                    features.append(feature)
                    file_names.append(filename)



                    #for i, feature_df in enumerate(features):

                    if evaluate == True:
                        #label function from train is used
                        feature_df_labeled = train.label_feature_df(feature, true_tract_data, discard_ambiguous=False, replicates=True)
                        #feature_df_labeled = train.label_feature_df_archie(feature_df, true_tracts_infer[i], discard_ambiguous=False, replicates=True)

                    else:

                        feature_df_labeled = label_feature_df_archie_infer_only_label(feature, true_tract_labeled, discard_ambiguous=False)


                    if not os.path.exists(os.path.join(output_dir, "feature_dfs")):
                        os.makedirs(os.path.join(output_dir, "feature_dfs"))

                    feature_df_labeled.to_csv(os.path.join(output_dir, "feature_dfs", "features" + str(filename) + ".csv") )

                    feature_df_labeleds.append(feature_df_labeled)


        replicate_counter = replicate_counter + 1


    
    #create one big test dataframe 
    test_df = pd.concat(feature_df_labeleds)

    #drop_dynamic_cols indicates whether all non-fixed size features should be dropped
    if drop_dynamic_cols == True:
        dynamic_cols = [col for col in test_df.columns if ('-ton' in col or col.startswith("pairwised_dist"))]
        test_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')

    #in the case of simulated data, we know the labels, but for inference we remove them for the test dataframe
    if simulated == True and average_for_inference == False:
        y_true = test_df["label"]
        test_df.drop(["label"], axis=1, inplace=True, errors='ignore')
    

    
    #statsmodel inference
    if statsmodel != None and compute_statsmodel==True:
        test_df_pred = test_df.copy()
        model_statsmodel = load_statsmodel(statsmodel)

        #only compute probabilities
        if compute_cutoffs == False:
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:       
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

            test_df_pred["probabilities"] = y_prob
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            new_dfs = split_dfs_sample_replicate(df_final)
            
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_wo_parallel(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            
            if plot_curves == True:
                plot_cutoffs(recalls, precisions)


        #statsmodel inference
    if statsmodel != None and compute_statsmodel==True:
        test_df_pred = test_df.copy()
        model_statsmodel = load_statsmodel(statsmodel)

        #only compute probabilities
        if compute_cutoffs == False:
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:       
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

            test_df_pred["probabilities"] = y_prob
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            new_dfs = split_dfs_sample_replicate(df_final)
            
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased_wo_parallel(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            
            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / Statsmodel unphased")
    

    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit)
        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_wo_parallel(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)

            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit phased")

            write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, "archie_1src_accuracy.txt")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            if evaluate == True:
                test_df_pred = test_df.copy()
                evaluate_scikit(test_df_pred, y_true, model_scikit)
       


    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit)
        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased_wo_parallel(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)

            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit unphased")

            write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, "archie_unphased_1src_accuracy.txt")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            #if evaluate == True:
            #    evaluate_scikit(test_df_pred, y_true, model_scikit)

    return test_df



global seq_len
global archaic_prop
global not_archaic_prop
global thread
global output_prefix
global output_dirs
global evaluate
global ref_ind_file
global tgt_ind_file
global anc_allele_file
global win_len
global win_step
global match_bonus
global max_mismatch 
global mismatch_penalty 
global process_archie

global discard_ambiguous
#def store_global_infer(seq_len, archaic_prop, not_archaic_prop, thread, output_prefix, output_dirs,  evaluate,  ref_ind_file, tgt_ind_file, anc_allele_file, win_len, win_step, match_bonus, max_mismatch, mismatch_penalty, process_archie):

from concurrent.futures import ProcessPoolExecutor as Pool
def store_global_infer( seq_len_new, archaic_prop_new, not_archaic_prop_new, thread_new, output_prefix_new, output_dirs_new,  evaluate_new,  ref_ind_file_new, tgt_ind_file_new, anc_allele_file_new, win_len_new, win_step_new, match_bonus_new, max_mismatch_new, mismatch_penalty_new, discard_ambiguous_new, process_archie_new):
    '''
    for using the pool for parallelization, the variables have to be stored globally
    '''
    global seq_len
    global archaic_prop
    global not_archaic_prop
    global thread
    global output_prefix
    global output_dirs
    global evaluate
    global ref_ind_file
    global tgt_ind_file
    global anc_allele_file
    global win_len
    global win_step
    global match_bonus
    global max_mismatch 
    global mismatch_penalty 
    global process_archie

    global discard_ambiguous

    seq_len = seq_len_new
    archaic_prop = archaic_prop_new
    not_archaic_prop = not_archaic_prop_new
    thread = thread_new
    output_prefix = output_prefix_new
    output_dirs = output_dirs_new
    evaluate = evaluate_new
    ref_ind_file = ref_ind_file_new
    tgt_ind_file = tgt_ind_file_new
    anc_allele_file = anc_allele_file_new
    win_len = win_len_new
    win_step = win_step_new
    match_bonus = match_bonus_new
    max_mismatch = max_mismatch_new
    mismatch_penalty = mismatch_penalty_new
    process_archie = process_archie_new

    discard_ambiguous = discard_ambiguous_new



def infer_parallel_process_label_data_sm(output_tuples):
    '''
    this function preprocesses the data for inference (i.e. calculating statistics)
    it is called via the pool parallelization process
    '''
    global seq_len
    global archaic_prop
    global not_archaic_prop
    global thread
    global output_prefix
    global output_dirs
    global evaluate
    global ref_ind_file
    global tgt_ind_file
    global anc_allele_file
    global win_len
    global win_step
    global match_bonus
    global max_mismatch 
    global mismatch_penalty 
    global process_archie


    file = output_tuples[1]
    replicate_counter = output_tuples[0]
    output_dir = output_tuples[2]

    
    #not necessary anymore
    #filename = os.path.splitext(file)[0]

    #true_tracts_infer = []
    feature_file = os.path.splitext(file)[0]+'.features'
    #computation of statistics
    preprocess.process_data(os.path.join(output_dir,file), ref_ind_file, tgt_ind_file, anc_allele_file, os.path.join(output_dir,feature_file), win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, process_archie)


    if file.endswith(".vcf"):
        true_tract = os.path.splitext(file)[0]+'.true.tracts.bed'
    else:
        true_tract = file.split(".")[0]+'.introgressed.tracts.bed'
    
    #reading not necessary
    #true_tract_data = pd.read_csv(os.path.join(output_dir, true_tract), sep="\t", header=None, names=['chr', 'start', 'end', 'hap', 'ind']) 
    #for infer necessary
    true_tract_data = pd.read_csv(os.path.join(output_dir, true_tract), sep="\t", header=None, names=['chr', 'start', 'end', 'hap', 'ind'])

    add_ind = None
    if true_tract_data is not None:
        true_tract_data["rep"] = replicate_counter
        if file.endswith(".vcf"):
            true_tract_data["hap"] = true_tract_data["hap"].str.replace("hap_", "")
        else:
        
            add_ind = "tsk_10"
            true_tract_data["hap"] = 1
            true_tract_data["ind"] = add_ind


    #labeling of true tracts
    true_tract_labeled = train._label(os.path.join(output_dir, true_tract), archaic_prop, not_archaic_prop, seq_len, add_ind)

    #label the true tracts according to the function from train
    if true_tract_labeled is not None:
        true_tract_labeled["rep"] = replicate_counter


    #load the feature files with statistics created before
    feature = pd.read_csv(os.path.join(output_dir, feature_file), sep="\t")
    feature["rep"] = replicate_counter


    if evaluate == True:
        #label function from train is used
        feature_df_labeled = train.label_feature_df(feature, true_tract_data, only_above_threshold=True, discard_ambiguous=False, replicates=True)
        #feature_df_labeled = train.label_feature_df_archie(feature_df, true_tracts_infer[i], discard_ambiguous=False, replicates=True)

    else:

        feature_df_labeled = label_feature_df_archie_infer_only_label(feature, true_tract_labeled, discard_ambiguous=False)


    if not os.path.exists(os.path.join(output_dir, "feature_dfs")):
        os.makedirs(os.path.join(output_dir, "feature_dfs"))

    return feature_df_labeled, true_tract_data


def predict_introgression_folders(nrep, nref, ntgt, seq_len, thread, output_prefix, output_dirs, statsmodel=None, scikitmodel=None, drop_dynamic_cols=True, evaluate=False, simulated=True, average_for_inference=False, compute_cutoffs=True, ref_ind_file=None, tgt_ind_file=None, model_name="archie", sample_name="sample1", compute_statsmodel = False, plot_curves=False, win_step_50k=False, discard_ambiguous=False, drop_kurtosis=False, drop_paired=False, drop_almost_all=False, use_haplotype_acc=False):    
    #set filenames for individuals, reference and target (needed for prprocess.process_data)

    if ref_ind_file == None:
        ref_ind_file = str(model_name) + "_new_infer" + "_nref" + str(nref) + "_ntgt" + str(ntgt) + ".ref.ind.list"
        tgt_ind_file  = str(model_name) + "_new_infer" + "_nref" + str(nref) + "_ntgt" + str(ntgt) + ".tgt.ind.list"
        train.create_ref_tgt_file(nref, ntgt, ref_ind_file, tgt_ind_file)


    anc_allele_file = None

    #set window length and stepsize
    
    win_len = 50000
    if win_step_50k == True:
        win_step = 50000
    else:
        win_step = 10000
 
    #I think these parameters are NOT necessary for ArchIE - just retained for the signature of preprocess.process_data
    match_bonus = 1
    max_mismatch = 1
    mismatch_penalty = 1

    process_archie = True

    #tracts with a proportion between not_archaic and archaic are labeled as ambiguous (in _label)
    archaic_prop = 0.7
    not_archaic_prop = 0.3
    

    true_tracts = []

    global true_tracts_infer
    true_tracts_infer = []

    features = []
    file_names = []
    replicate_counter = 0

    feature_df_labeleds = []

    #new parallel part
    store_global_infer(seq_len, archaic_prop, not_archaic_prop, thread, output_prefix, output_dirs,  evaluate,  ref_ind_file, tgt_ind_file, anc_allele_file, win_len, win_step, match_bonus, max_mismatch, mismatch_penalty, discard_ambiguous, process_archie)


    output_tuples = []

    replicate_counter = 0
    for output_dir in output_dirs:
        if os.path.isdir(output_dir):
            for replicate, file in enumerate(os.listdir(output_dir)):
                if file.endswith(".vcf") or file.endswith("biallelic.vcf.gz"):
                    output_tuples.append((replicate_counter, file, output_dir))
                    replicate_counter = replicate_counter + 1


    feature_df_labeleds = []
    pool = Pool()      

    feature_df_labeleds, true_tracts_infer = zip(* pool.map(infer_parallel_process_label_data_sm, output_tuples) )

    #create one big test dataframe 
    test_df = pd.concat(feature_df_labeleds)

    true_tracts_infer =list(true_tracts_infer)

    #drop_dynamic_cols indicates whether all non-fixed size features should be dropped
    if drop_dynamic_cols == True:
        dynamic_cols = [col for col in test_df.columns if ('-ton' in col or col.startswith("pairwised_dist"))]
        test_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')


    #also for the other 'reduced forms'

    if drop_kurtosis == True:
        dynamic_cols = [col for col in test_df.columns if ('kurtosis_pairwised_dist' in col or col.startswith("pairwised_dist")) ]
        test_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')

    if drop_paired == True:
        dynamic_cols =  [col for col in test_df.columns if (col.startswith("pairwised_dist"))]
        test_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')

    if drop_almost_all == True:
        dynamic_cols = [col for col in test_df.columns if ('-ton' in col or 'pairwised_dist' in col )]
        test_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')

    #in the case of simulated data, we know the labels, but for inference we remove them for the test dataframe
    if simulated == True and average_for_inference == False:
        y_true = test_df["label"]
        test_df.drop(["label"], axis=1, inplace=True, errors='ignore')
    
    global new_dfs
    
    #statsmodel inference
    if statsmodel != None and compute_statsmodel==True:
        test_df_pred = test_df.copy()
        model_statsmodel = load_statsmodel(statsmodel)

        #only compute probabilities
        if compute_cutoffs == False:
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:       
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

            test_df_pred["probabilities"] = y_prob
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            new_dfs = split_dfs_sample_replicate(df_final)
            
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            
            if plot_curves == True:
                plot_cutoffs(recalls, precisions)


        #statsmodel inference
    if statsmodel != None and compute_statsmodel==True:
        test_df_pred = test_df.copy()
        model_statsmodel = load_statsmodel(statsmodel)

        #only compute probabilities
        if compute_cutoffs == False:
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:       
            y_prob = infer_statsmodel(test_df_pred.copy(), model_statsmodel)

            test_df_pred["probabilities"] = y_prob
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            new_dfs = split_dfs_sample_replicate(df_final)
            
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            
            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / Statsmodel unphased")
    

    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            pass
            '''
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.txt", model_name + " ,rep: " + str(nrep) + " ,nref: " + str(nref) + " ,ntgt: " + str(ntgt))
            '''

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
           
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]

            #new parallel
            pool = Pool()      

            if use_haplotype_acc == False:
                precisions, recalls = zip(* pool.map(cal_accuracy_v2, cut_offs) )
            else:
                precisions, recalls = zip(* pool.map(cal_accuracy_haplotypes_phased, cut_offs) )

            '''
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            '''

            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit phased")


            if use_haplotype_acc == False:
                write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, model_name +  sample_name + "archie_1src_accuracy.txt")
            else:
                write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, model_name +  sample_name + "archie_1src_accuracy_haplotypes.txt")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            '''
            if evaluate == True:
                test_df_pred = test_df.copy()
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.txt", model_name + " ,rep: " + str(nrep) + " ,nref: " + str(nref) + " ,ntgt: " + str(ntgt))
            '''


    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.png", model_name + ", rep: " + str(nrep) + ", nref: " + str(nref) + ", ntgt: " + str(ntgt))
        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]

            #new parallel
            pool = Pool()      
            

            if use_haplotype_acc == False:
                precisions, recalls = zip(* pool.map(cal_accuracy_v2_unphased, cut_offs) )
            else:
                precisions, recalls = zip(* pool.map(cal_accuracy_haplotypes_unphased, cut_offs) )


            '''
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            '''
            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit unphased")

            if use_haplotype_acc == False:
                write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, model_name +  sample_name + "archie_unphased_1src_accuracy.txt")
            else:
                write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, model_name +  sample_name + "archie_unphased_1src_accuracy_haplotypes.txt")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            if evaluate == True:
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.png", model_name + " ,rep: " + str(nrep) + " ,nref: " + str(nref) + " ,ntgt: " + str(ntgt))
    return test_df





def predict_introgression_folders_allmodeltypes(nrep, nref, ntgt, seq_len, thread, output_prefix, output_dirs, statsmodel=None, scikitmodel=None, drop_dynamic_cols=True, evaluate=False, simulated=True, average_for_inference=False, compute_cutoffs=True, ref_ind_file=None, tgt_ind_file=None, model_name="archie", sample_name="sample1", compute_statsmodel = False, plot_curves=False, win_step_50k=False, discard_ambiguous=False, drop_kurtosis=True, drop_paired=True, drop_almost_all=True):    
    #set filenames for individuals, reference and target (needed for prprocess.process_data)

    if ref_ind_file == None:
        ref_ind_file = str(model_name) + "_new_infer" + "_nref" + str(nref) + "_ntgt" + str(ntgt) + ".ref.ind.list"
        tgt_ind_file  = str(model_name) + "_new_infer" + "_nref" + str(nref) + "_ntgt" + str(ntgt) + ".tgt.ind.list"
        train.create_ref_tgt_file(nref, ntgt, ref_ind_file, tgt_ind_file)


    anc_allele_file = None

    #set window length and stepsize
    
    win_len = 50000
    if win_step_50k == True:
        win_step = 50000
    else:
        win_step = 10000
 
    #I think these parameters are NOT necessary for ArchIE - just retained for the signature of preprocess.process_data
    match_bonus = 1
    max_mismatch = 1
    mismatch_penalty = 1

    process_archie = True

    #tracts with a proportion between not_archaic and archaic are labeled as ambiguous (in _label)
    archaic_prop = 0.7
    not_archaic_prop = 0.3
    

    true_tracts = []

    global true_tracts_infer
    true_tracts_infer = []

    features = []
    file_names = []
    replicate_counter = 0

    feature_df_labeleds = []

    #new parallel part
    store_global_infer(seq_len, archaic_prop, not_archaic_prop, thread, output_prefix, output_dirs,  evaluate,  ref_ind_file, tgt_ind_file, anc_allele_file, win_len, win_step, match_bonus, max_mismatch, mismatch_penalty, discard_ambiguous, process_archie)


    output_tuples = []

    replicate_counter = 0
    for output_dir in output_dirs:
        if os.path.isdir(output_dir):
            for replicate, file in enumerate(os.listdir(output_dir)):
                if file.endswith(".vcf") or file.endswith("biallelic.vcf.gz"):
                    output_tuples.append((replicate_counter, file, output_dir))
                    replicate_counter = replicate_counter + 1


    feature_df_labeleds = []
    pool = Pool()      

    feature_df_labeleds, true_tracts_infer = zip(* pool.map(infer_parallel_process_label_data_sm, output_tuples) )

    #create one big test dataframe 
    test_df = pd.concat(feature_df_labeleds)

    true_tracts_infer =list(true_tracts_infer)


    #in the case of simulated data, we know the labels, but for inference we remove them for the test dataframe
    if simulated == True and average_for_inference == False:
        y_true = test_df["label"]
        test_df.drop(["label"], axis=1, inplace=True, errors='ignore')


    scikit_full_inference(test_df, y_true, nrep, nref, ntgt, scikitmodel, compute_cutoffs, plot_curves, model_name, sample_name, type_name="")


    #drop_dynamic_cols indicates whether all non-fixed size features should be dropped
    if drop_dynamic_cols == True:
        dynamic_cols = [col for col in test_df.columns if ('-ton' in col or col.startswith("pairwised_dist"))]
        test_df_drop_dynamic_cols = test_df.copy()
        test_df_drop_dynamic_cols.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')
        scikit_full_inference(test_df, y_true, nrep, nref, ntgt, "fixed_" + scikitmodel, compute_cutoffs, plot_curves, model_name, sample_name, type_name="_dropdyn_")


    #also for the other 'reduced forms'

    if drop_kurtosis == True:
        dynamic_cols = [col for col in test_df.columns if ('kurtosis_pairwised_dist' in col or col.startswith("pairwised_dist")) ]
        test_df_drop_dynamic_cols = test_df.copy()
        test_df_drop_dynamic_cols.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')
        scikit_full_inference(test_df, y_true, nrep, nref, ntgt, "nokurt_" + scikitmodel, compute_cutoffs, plot_curves, model_name, sample_name, type_name="_dropkurt_")

    if drop_paired == True:
        dynamic_cols =  [col for col in test_df.columns if (col.startswith("pairwised_dist"))]
        test_df_drop_dynamic_cols = test_df.copy()
        test_df_drop_dynamic_cols.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')
        scikit_full_inference(test_df, y_true, nrep, nref, ntgt, "nopaired_" + scikitmodel, compute_cutoffs, plot_curves, model_name, sample_name, type_name="_droppaired_")

    if drop_almost_all == True:
        dynamic_cols = [col for col in test_df.columns if ('-ton' in col or 'pairwised_dist' in col )]
        test_df_drop_dynamic_cols = test_df.copy()
        test_df_drop_dynamic_cols.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')
        scikit_full_inference(test_df, y_true, nrep, nref, ntgt, "fullreduced_" + scikitmodel, compute_cutoffs, plot_curves, model_name, sample_name, type_name="_dropalmostall_")

    


def scikit_full_inference(test_df, y_true, nrep, nref, ntgt, scikitmodel, compute_cutoffs, plot_curves, model_name, sample_name, type_name=""):
    '''
    this function creates precision-recall-curves using scikit
    '''
    global true_tracts_infer
    global new_dfs

    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            pass
            '''
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.txt", model_name + " ,rep: " + str(nrep) + " ,nref: " + str(nref) + " ,ntgt: " + str(ntgt))
            '''

        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
           
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]

            #new parallel
            pool = Pool()      


            precisions, recalls = zip(* pool.map(cal_accuracy_v2, cut_offs) )

            '''
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            '''

            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit phased")

            write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, model_name +  sample_name + type_name + "archie_1src_accuracy.txt")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            '''
            if evaluate == True:
                test_df_pred = test_df.copy()
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.txt", model_name + " ,rep: " + str(nrep) + " ,nref: " + str(nref) + " ,ntgt: " + str(ntgt))
            '''


    #scikit inference
    if scikitmodel != None:
        test_df_pred = test_df.copy()
        model_scikit = load_scikit(scikitmodel)

        if compute_cutoffs == False:
            y_prob = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
            if evaluate == True:
                #create precision-recall curve using all instances from dataframe
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + type_name + "archie_scikit_precrec.png", model_name + ", rep: " + str(nrep) + ", nref: " + str(nref) + ", ntgt: " + str(ntgt))
        #create cut-offs and compute precision-recall curve as in sstar-analysis
        else:
            test_df_pred.drop(["probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
            y_pred = infer_scikit(test_df_pred.copy(), model_scikit, probabilities=True)
              
            test_df_pred["probabilities"] = y_pred
            #new column, unneccessary should be changed
            #start and end are necessary for comparing true and inferred tracts (see cal_accuracy_v2)
            test_df_pred["tractstart"] = test_df_pred["start"]
            test_df_pred["tractend"] = test_df_pred["end"]

            df_final = test_df_pred
            
            #create a list of dataframes, each dataframe corresponds to one individual / one haplotype / one replicate
            new_dfs = split_dfs_sample_replicate(df_final)
           
            precisions = []
            recalls = []
            cut_offs = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999]

            #new parallel
            pool = Pool()      


            precisions, recalls = zip(* pool.map(cal_accuracy_v2_unphased, cut_offs) )


            '''
            for cut_off in cut_offs:
                prec, rec = cal_accuracy_v2_unphased(true_tracts_infer, new_dfs, cutoff = cut_off)
                precisions.append(prec)
                recalls.append(rec)
            '''
            if plot_curves == True:
                plot_cutoffs(recalls, precisions, title="Precision-Recall curve for computed cutoffs / scikit unphased")

            write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, model_name +  sample_name + "archie_unphased_1src_accuracy.txt")
   
            #scikit_precision_recall(df_final, plot_label="scikitmodel")

            if evaluate == True:
                evaluate_scikit(test_df_pred, y_true, model_scikit, plot_curves, True, model_name +  sample_name + "archie_scikit_precrec.png", model_name + " ,rep: " + str(nrep) + " ,nref: " + str(nref) + " ,ntgt: " + str(ntgt))
    return test_df






def write_prec_recall_df(cut_offs, precisions, recalls, model_name, sample_name, acc_filename):
    '''
    write precision-recall-information to file
    '''
    
    prec_rec_list = []
    for i, cut_off in enumerate(cut_offs):
        prec_rec_list.append([model_name, sample_name, cut_off, precisions[i], recalls[i]])

    prec_rec_df = pd.DataFrame(prec_rec_list)
    prec_rec_df.columns = ["demography","sample","cutoff","precision","recall"]

    if not os.path.exists(os.path.join("results", "inference", "archie")):
            os.makedirs(os.path.join("results", "inference", "archie"))


    prec_rec_df.to_csv(os.path.join("results", "inference", "archie",acc_filename), sep="\t", index=False, na_rep="nan")
    


def plot_cutoffs(recs, precs, title=None):
    '''
    simple plot function for precision-recall
    '''
    plt.plot(recs, precs)
    plt.scatter(recs, precs)
    plt.xlim(left=0)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision-Recall curve for computed cutoffs")
    if title != None:
        plt.title(title)
    plt.show()


def create_testdata(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed=None):
    '''
    create testdata for inference - one folder with many files is created
    '''
    train._simulation_manager(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed)


def create_testdata_folders(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed=None):
    '''
    create testdata for inference - one folder for each replicate is created
    '''
    train._simulation_manager_folders(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed)

def compute_average_probabilities_true_tracts(df, probabilities, true_tracts, output_file=None, cut_off=0.5):
    #NOT used at the moment
    #the idea was to compute the probabilities over one overlap region, probably not necessary

    df["probabilities"] = probabilities

    if output_file is not None:
        df.to_csv(output_file)

    inferred_tracts = []

    for tract_sample in true_tracts:
        for entry in tract_sample.values.tolist():

            start = entry[1]
            end = entry[2]

            haplo = int(entry[3])
            ind = entry[4]

            replicate = entry[5]

            inferred_df = df.copy()
         
            inferred_df = inferred_df[(inferred_df["sample"] == ind) & (inferred_df["haplo"] == haplo) & (inferred_df["rep"] == replicate)]
          

            inferred_df.drop(["overlap"], axis=1, inplace=True, errors='ignore')

            inferred_df['overlap2'] = 0
            inferred_df["overlap2"] = inferred_df.apply(lambda row: row.overlap2 + train.getOverlap(int(start),int(end), int(row.start),int(row.end)) , axis=1)

            inferred_df = inferred_df[inferred_df["overlap2"] > 0]

            #average over probabilities
            inferred_df = inferred_df.groupby(by=['chrom','haplo', 'sample', 'rep'])['probabilities'].mean().reset_index()

        
            inferred_df["tractstart"] = int(start)
            inferred_df["tractend"] = int(end)

            inferred_tracts.append(inferred_df.copy())

    inferred_tracts = pd.concat(inferred_tracts)

    inferred_tracts.drop(['overlap_length', 'interval', 'overlap_percentage', "label_one_1", "label_one_2", "label_one_3"], axis=1, inplace=True, errors='ignore')
    return inferred_tracts



def load_statsmodel(model_file):
    '''
    load statsmodel model
    '''
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    model = sm.load(model_file)

    return model


def load_scikit(model_file):
    '''
    load scikit model
    '''
    import pickle
    with open(model_file , 'rb') as f:
        model_scikit = pickle.load(f)
    return model_scikit


def infer_statsmodel(X_test, model):
    '''
    compute probabilities for test set with statsmodels logistic classification
    '''
    import statsmodels.api as sm
    X_test.replace(np.nan, 0, inplace=True)
    X_test.replace(pd.NA, 0, inplace=True)

    X_test.drop(['overlap_length','chrom', 'start', 'end', 'interval', 'overlap_percentage', "label_one_1", "label_one_2", "label_one_3"], axis=1, inplace=True, errors='ignore')

    X_test.drop(['overlap','sample', 'haplo', 'rep'], axis=1, inplace=True, errors='ignore')

    X_test.drop(["label"], axis=1, inplace=True, errors='ignore')


    X_test = X_test.astype(float)


    predictions = model.predict(sm.add_constant(X_test, prepend=False))

    return predictions
    

def infer_scikit(X_test, model, evaluate = True, probabilities = False):
    '''
    compute probabilities for test set with scikit logistic classification
    '''
    X_test.replace(np.nan, 0, inplace=True)
    X_test.replace(pd.NA, 0, inplace=True)
    
    X_test.drop(['overlap_length','chrom', 'start', 'end', 'interval', 'overlap_percentage', "label_one_1", "label_one_2", "label_one_3"], axis=1, inplace=True, errors='ignore')
    X_test.drop(['overlap','sample', 'haplo', 'rep'], axis=1, inplace=True, errors='ignore')
    X_test.drop(["label"], axis=1, inplace=True, errors='ignore')

    if probabilities == False:
        y_pred = model.predict(X_test)
    else:
        y_pred = model.predict_proba(X_test)

    return y_pred[:,1]


def scikit_precision_recall(X_test, plot_label="prediction logistic/scikit"):
    '''
    given a dataframe with labels and predicted probabilities, create a precision-recall curve via scikit
    '''
    from sklearn.metrics import precision_recall_curve

    X_test.drop(['overlap_length','chrom', 'start', 'end', 'interval', 'overlap_percentage', "label_one_1", "label_one_2", "label_one_3"], axis=1, inplace=True, errors='ignore')


    y_true = X_test["label"]
    X_test.drop(["label"], axis=1, inplace=True, errors='ignore')

    lr_probs = X_test["probabilities"]
    X_test.drop(["probabilities"], axis=1, inplace=True, errors='ignore')

    precision, recall, thresholds = precision_recall_curve(y_true.astype(int), lr_probs)

    plt.plot(recall, precision, marker='.', label=plot_label)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.show()


def evaluate_scikit(X_test, y, model, plot_curves = False, write_model=True, filename=None, textlabel = "prediction logistic/scikit" ):
    '''
    given a dataframe and a list of the according labels, predict probabilities for the instances of the df and plot a precision-recall curve
    '''
    from sklearn.metrics import precision_recall_curve
    
    X_test.drop(['haplo', 'rep', 'sample', 'overlap_length','chrom', 'start', 'end', 'interval', 'overlap_percentage', "label_one_1", "label_one_2", "label_one_3"], axis=1, inplace=True, errors='ignore')
    X_test.drop(["overlap", "probabilities", "tractstart", "tractend"], axis=1, inplace=True, errors='ignore')
    X_test.drop(["label"], axis=1, inplace=True, errors='ignore')
    

    lr_probs = model.predict_proba(X_test)

    print(np.unique(y))

    precision, recall, thresholds = precision_recall_curve(y.astype(int), lr_probs[:,1])

    new_precrec_df = pd.DataFrame([precision, recall, thresholds])
    
    
    plt.plot(recall, precision, marker='.', label=textlabel)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()

    if write_model == True and filename != None:
        if not os.path.exists(os.path.join("results", "inference", "archie")):
            os.makedirs(os.path.join("results", "inference", "archie"))
        plt.savefig(os.path.join("results", "inference", "archie",filename))
        #new_precrec_df.to_csv(os.path.join("results", "inference", "archie",filename), sep="\t", index=False, na_rep="nan")

    if plot_curves == True:
        plt.show()





def label_feature_df_archie_infer_only_label(feature_df, true_tract_list, discard_ambiguous=False):
    """
    Description:
        Helper function to label a fragment as 'introgressed' or 'not introgressed'
        If one evaluates the true tracts in the calculation of precision and recall, this is no necessary;
        however, using discard_ambiguous, one can remove all ambiguous windows

    Arguments:
        feature_df DataFrame: Dataframe containing all windows/fragments to be labeled
        true_tract_list DataFrame: Dataframe containing all true tracts which are used for labeling
        discard_ambiguous boolean: discard windows/fragments of intermediate introgression content (according to true_tract_list)
    """

    feature_df['label'] = 0

    if true_tract_list is None:
        return feature_df

    
    true_tract_list["hap"] = true_tract_list["hap"].astype(str).str.replace("hap_", "")
    true_tract_list["hap"] = true_tract_list["hap"].astype(int)

    for ie, entry in enumerate(true_tract_list.values.tolist()):
   
        haplo = int(entry[0])
        ind = entry[1]
        tract_label = entry[4]
        replicate = entry[5]
      
        conditions = (feature_df["sample"] == ind) & (feature_df["haplo"] == haplo) & (feature_df["rep"] == replicate)
     

        if tract_label[0] == 1:
            feature_df.loc[conditions, "label"] = 1

        
        if tract_label[2] == 1:
            if discard_ambiguous == True:
            
                 feature_df = feature_df.drop(feature_df[conditions].index)
     
    return feature_df



def cal_accuracy_v2(cutoff=0.5):
    import pybedtools
    import numpy as np
    import pandas as pd
    from copy import deepcopy
    """
    Description:
        Helper function for calculating accuracy; in contrast to cal_accuracy from utils, it is iterated over all replicates and samples.

    Arguments:
        true_tracts list: list of dataframes containing true introgresssed tracts.
        inferred_dfs DataFrame: dataframe containing inferred tracts

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """

    true_tracts = pd.concat(deepcopy(true_tracts_infer))

    all_total_inferred_tracts = 0
    all_total_true_tracts = 0
    all_true_positives = 0

    true_tracts["processed"] = 0

    for replicate in new_dfs:
            
        for sample in replicate:

            for haplo in sample:
                haplo = pd.DataFrame(haplo)

                #apply cutoff
                haplo = haplo[haplo["probabilities"] >= cutoff]

                if haplo.shape[0] > 0:
                
                    curr_rep = int(haplo["rep"].unique()[0])
                    curr_sample = haplo["sample"].unique()[0]
                    curr_haplo = int(haplo["haplo"].unique()[0])
                    
                    conditions = (true_tracts["ind"].astype(str) == curr_sample) & (true_tracts["hap"].astype(int) == curr_haplo) & (true_tracts["rep"].astype(int) == curr_rep)
                    curr_truth_tracts = true_tracts[conditions]

                    true_tracts.loc[conditions, "processed"] = 1
                
                    truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
                
                    inferred_tracts = haplo[["chrom", "tractstart", "tractend"]]
                    inferred_tracts.columns = ["chrom", "start", "end"]
                    truth_tracts.columns = ["chrom", "start", "end"]

                    truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
                    inferred_tracts = pybedtools.BedTool.from_dataframe(inferred_tracts).sort().merge()

                    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
                    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
                    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

                    all_total_inferred_tracts = all_total_inferred_tracts + total_inferred_tracts
                    all_total_true_tracts = all_total_true_tracts + total_true_tracts
                    all_true_positives = all_true_positives + true_positives
                        
    true_tracts = true_tracts[true_tracts["processed"] == 0]
    for i, tract in true_tracts.iterrows():
        curr_truth_tracts = tract.to_frame().T
        truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
        truth_tracts.columns = ["chrom", "start", "end"]
        truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        all_total_true_tracts = all_total_true_tracts + total_true_tracts
        if cutoff == 0:
            all_true_positives = all_true_positives + total_true_tracts



    if float(all_total_inferred_tracts) == 0: precision = np.nan
    else: precision = all_true_positives / float(all_total_inferred_tracts) * 100
    if float(all_total_true_tracts) == 0: recall = np.nan
    else: recall = all_true_positives / float(all_total_true_tracts) * 100

    return precision, recall


def cal_accuracy_v2_unphased(cutoff=0.5):
    import pybedtools
    import numpy as np
    import pandas as pd
    from copy import deepcopy
    """
    Description:
        Helper function for calculating accuracy; in contrast to cal_accuracy from utils, it is iterated over all replicates and samples.

    Arguments:
        true_tracts list: list of dataframes containing true introgresssed tracts.
        inferred_dfs DataFrame: dataframe containing inferred tracts

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """

    true_tracts = pd.concat(deepcopy(true_tracts_infer))

    all_total_inferred_tracts = 0
    all_total_true_tracts = 0
    all_true_positives = 0

    true_tracts["processed"] = 0

    for replicate in new_dfs:
            
        for sample in replicate:
            haplos = []
            for haplo in sample:
                haplos.append(pd.DataFrame(haplo))
            haplo = pd.concat(haplos)

            haplo = pd.DataFrame(haplo)

            #apply cutoff
            haplo = haplo[haplo["probabilities"] >= cutoff]

            if haplo.shape[0] > 0:
            
                curr_rep = int(haplo["rep"].unique()[0])
                curr_sample = haplo["sample"].unique()[0]
                
                conditions = (true_tracts["ind"].astype(str) == curr_sample) & (true_tracts["rep"].astype(int) == curr_rep)
                curr_truth_tracts = true_tracts[conditions]

                true_tracts.loc[conditions, "processed"] = 1
            
                truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
            
                inferred_tracts = haplo[["chrom", "tractstart", "tractend"]]
                inferred_tracts.columns = ["chrom", "start", "end"]
                truth_tracts.columns = ["chrom", "start", "end"]

                truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
                inferred_tracts = pybedtools.BedTool.from_dataframe(inferred_tracts).sort().merge()

                total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
                total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
                true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

                all_total_inferred_tracts = all_total_inferred_tracts + total_inferred_tracts
                all_total_true_tracts = all_total_true_tracts + total_true_tracts
                all_true_positives = all_true_positives + true_positives
                        
    true_tracts = true_tracts[true_tracts["processed"] == 0]

    haplo_true_tracts = true_tracts.groupby(["rep", "ind"])

    for i, haplo_group in haplo_true_tracts:

        truth_tracts = haplo_group[["chr", "start", "end"]]
        truth_tracts.columns = ["chrom", "start", "end"]
        truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        all_total_true_tracts = all_total_true_tracts + total_true_tracts
        if cutoff == 0:
            all_true_positives = all_true_positives + total_true_tracts



    if float(all_total_inferred_tracts) == 0: precision = np.nan
    else: precision = all_true_positives / float(all_total_inferred_tracts) * 100
    if float(all_total_true_tracts) == 0: recall = np.nan
    else: recall = all_true_positives / float(all_total_true_tracts) * 100

    return precision, recall



def cal_accuracy_haplotypes_unphased(cutoff=0.5, seqlen = 50000, archaic_prop = 0.7, not_archaic_prop = 0.3):
    import pybedtools
    import numpy as np
    import pandas as pd
    from copy import deepcopy
    """
    Description:
        Helper function for calculating accuracy; in contrast to cal_accuracy from utils, it is iterated over all replicates and samples.

    Arguments:
        true_tracts list: list of dataframes containing true introgresssed tracts.
        inferred_dfs DataFrame: dataframe containing inferred tracts

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """

    true_tracts = pd.concat(deepcopy(true_tracts_infer))


    true_label_list = []
    infer_label_list = []

    true_tracts["processed"] = 0

    for replicate in new_dfs:
            
        for sample in replicate:
            haplos = []
            for haplo in sample:
                haplos.append(pd.DataFrame(haplo))
            haplo = pd.concat(haplos)

            haplo = pd.DataFrame(haplo)

            #apply cutoff
            haplo = haplo[haplo["probabilities"] >= cutoff]

            if haplo.shape[0] > 0:
            
                curr_rep = int(haplo["rep"].unique()[0])
                curr_sample = haplo["sample"].unique()[0]
                
                conditions = (true_tracts["ind"].astype(str) == curr_sample) & (true_tracts["rep"].astype(int) == curr_rep)
                curr_truth_tracts = true_tracts[conditions]

                true_tracts.loc[conditions, "processed"] = 1
            
                truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
            
                inferred_tracts = haplo[["chrom", "tractstart", "tractend"]]
                inferred_tracts.columns = ["chrom", "start", "end"]
                truth_tracts.columns = ["chrom", "start", "end"]

                truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
                inferred_tracts = pybedtools.BedTool.from_dataframe(inferred_tracts).sort().merge()

                total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
                total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
                true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])


                if (total_true_tracts <= int(not_archaic_prop * seqlen)) or (total_true_tracts >= int(archaic_prop * seqlen)):

            
                    if total_true_tracts <= int(not_archaic_prop * seqlen):
                        true_label_list.append(0)

                    if total_true_tracts >= int(archaic_prop * seqlen):
                        true_label_list.append(1) 

                    if total_inferred_tracts <= int(not_archaic_prop * seqlen):
                        infer_label_list.append(0)

                    elif total_inferred_tracts >= int(archaic_prop * seqlen):
                        infer_label_list.append(1) 

                    else:
                        infer_label_list.append(2)


                        
    true_tracts = true_tracts[true_tracts["processed"] == 0]

    haplo_true_tracts = true_tracts.groupby(["rep", "ind"])

    for i, haplo_group in haplo_true_tracts:
        truth_tracts = haplo_group[["chr", "start", "end"]]

        truth_tracts.columns = ["chrom", "start", "end"]
        truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()

        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        if (total_true_tracts <= int(not_archaic_prop * seqlen)) or (total_true_tracts >= int(archaic_prop * seqlen)):
            if total_true_tracts <= int(not_archaic_prop * seqlen):
                true_label_list.append(0)
                infer_label_list.append(0)

            if total_true_tracts >= int(archaic_prop * seqlen):
                true_label_list.append(1) 
                infer_label_list.append(0)



    true_label_array = np.array(true_label_list)
    infer_label_array = np.array(infer_label_list)

    from sklearn.metrics import recall_score
    from sklearn.metrics import precision_score

    precision = precision_score(true_label_array, infer_label_array)
    recall = recall_score(true_label_array, infer_label_array)


    return precision, recall


def cal_accuracy_haplotypes_phased(cutoff=0.5, seqlen = 50000, archaic_prop = 0.7, not_archaic_prop = 0.3):
    import pybedtools
    import numpy as np
    import pandas as pd
    from copy import deepcopy
    """
    Description:
        Helper function for calculating accuracy; in contrast to cal_accuracy from utils, it is iterated over all replicates and samples.

    Arguments:
        true_tracts list: list of dataframes containing true introgresssed tracts.
        inferred_dfs DataFrame: dataframe containing inferred tracts

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """

    true_tracts = pd.concat(deepcopy(true_tracts_infer))


    true_label_list = []
    infer_label_list = []

    true_tracts["processed"] = 0

    for replicate in new_dfs:
            
        for sample in replicate:

            for haplo in sample:
                haplo = pd.DataFrame(haplo)

                #apply cutoff
                haplo = haplo[haplo["probabilities"] >= cutoff]

                if haplo.shape[0] > 0:

            
                    curr_rep = int(haplo["rep"].unique()[0])
                    curr_sample = haplo["sample"].unique()[0]
                    curr_haplo = int(haplo["haplo"].unique()[0])
                    
                    conditions = (true_tracts["ind"].astype(str) == curr_sample) & (true_tracts["hap"].astype(int) == curr_haplo)  & (true_tracts["rep"].astype(int) == curr_rep)
                    curr_truth_tracts = true_tracts[conditions]

                    true_tracts.loc[conditions, "processed"] = 1
                
                    truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
                
                    inferred_tracts = haplo[["chrom", "tractstart", "tractend"]]
                    inferred_tracts.columns = ["chrom", "start", "end"]
                    truth_tracts.columns = ["chrom", "start", "end"]

                    truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
                    inferred_tracts = pybedtools.BedTool.from_dataframe(inferred_tracts).sort().merge()

                    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
                    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
                    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])


                    if (total_true_tracts <= int(not_archaic_prop * seqlen)) or (total_true_tracts >= int(archaic_prop * seqlen)):

                
                        if total_true_tracts <= int(not_archaic_prop * seqlen):
                            true_label_list.append(0)

                        if total_true_tracts >= int(archaic_prop * seqlen):
                            true_label_list.append(1) 

                        if total_inferred_tracts <= int(not_archaic_prop * seqlen):
                            infer_label_list.append(0)

                        elif total_inferred_tracts >= int(archaic_prop * seqlen):
                            infer_label_list.append(1) 

                        else:
                            infer_label_list.append(2)

                        
    true_tracts = true_tracts[true_tracts["processed"] == 0]
    for i, tract in true_tracts.iterrows():
        curr_truth_tracts = tract.to_frame().T
        truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
        truth_tracts.columns = ["chrom", "start", "end"]
        truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()

        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        if (total_true_tracts <= int(not_archaic_prop * seqlen)) or (total_true_tracts >= int(archaic_prop * seqlen)):
            if total_true_tracts <= int(not_archaic_prop * seqlen):
                true_label_list.append(0)
                infer_label_list.append(0)

            if total_true_tracts >= int(archaic_prop * seqlen):
                true_label_list.append(1) 
                infer_label_list.append(0)



    true_label_array = np.array(true_label_list)
    infer_label_array = np.array(infer_label_list)

    from sklearn.metrics import recall_score
    from sklearn.metrics import precision_score

    precision = precision_score(true_label_array, infer_label_array)
    recall = recall_score(true_label_array, infer_label_array)


    return precision, recall


def cal_accuracy_v2_wo_parallel(true_tracts, inferred_dfs, cutoff=0.5):
    import pybedtools
    import numpy as np
    import pandas as pd
    from copy import deepcopy
    """
    Description:
        Helper function for calculating accuracy; in contrast to cal_accuracy from utils, it is iterated over all replicates and samples.
        This version is not used anymore 

    Arguments:
        true_tracts list: list of dataframes containing true introgresssed tracts.
        inferred_dfs DataFrame: dataframe containing inferred tracts

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """

    true_tracts = pd.concat(deepcopy(true_tracts))

    all_total_inferred_tracts = 0
    all_total_true_tracts = 0
    all_true_positives = 0


    for replicate in inferred_dfs:
            
        for sample in replicate:

            for haplo in sample:
                haplo = pd.DataFrame(haplo)

                #apply cutoff
                haplo = haplo[haplo["probabilities"] >= cutoff]

                if haplo.shape[0] > 0:
                
                    curr_rep = int(haplo["rep"].unique()[0])
                    curr_sample = haplo["sample"].unique()[0]
                    curr_haplo = int(haplo["haplo"].unique()[0])
                    
                    conditions = (true_tracts["ind"].astype(str) == curr_sample) & (true_tracts["hap"].astype(int) == curr_haplo) & (true_tracts["rep"].astype(int) == curr_rep)
                    curr_truth_tracts = true_tracts[conditions]
                    true_tracts.drop(true_tracts[conditions].index)
                
                    truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
                
                    inferred_tracts = haplo[["chrom", "tractstart", "tractend"]]
                    inferred_tracts.columns = ["chrom", "start", "end"]
                    truth_tracts.columns = ["chrom", "start", "end"]

                    truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
                    inferred_tracts = pybedtools.BedTool.from_dataframe(inferred_tracts).sort().merge()

                    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
                    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
                    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

                    all_total_inferred_tracts = all_total_inferred_tracts + total_inferred_tracts
                    all_total_true_tracts = all_total_true_tracts + total_true_tracts
                    all_true_positives = all_true_positives + true_positives
                        
    
    for i, tract in true_tracts.iterrows():
        curr_truth_tracts = tract.to_frame().T
        truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
        truth_tracts.columns = ["chrom", "start", "end"]
        truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        all_total_true_tracts = all_total_true_tracts + total_true_tracts


    if float(all_total_inferred_tracts) == 0: precision = np.nan
    else: precision = all_true_positives / float(all_total_inferred_tracts) * 100
    if float(all_total_true_tracts) == 0: recall = np.nan
    else: recall = all_true_positives / float(all_total_true_tracts) * 100

    return precision, recall



def cal_accuracy_v2_unphased_wo_parallel(true_tracts, inferred_dfs, cutoff=0.5):
    import pybedtools
    import numpy as np
    import pandas as pd
    from copy import deepcopy
    """
    Description:
        Helper function for calculating accuracy; in contrast to cal_accuracy from utils, it is iterated over all replicates and samples.
        This version is not used anymore 

    Arguments:
        true_tracts list: list of dataframes containing true introgresssed tracts.
        inferred_dfs DataFrame: dataframe containing inferred tracts

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """

    true_tracts = pd.concat(deepcopy(true_tracts))

    all_total_inferred_tracts = 0
    all_total_true_tracts = 0
    all_true_positives = 0

    for replicate in inferred_dfs:
            
        for sample in replicate:
            haplos = []
            for haplo in sample:
                haplos.append(pd.DataFrame(haplo))
            haplo = pd.concat(haplos)

            #apply cutoff
            haplo = haplo[haplo["probabilities"] >= cutoff]

            if haplo.shape[0] > 0:
            
                curr_rep = int(haplo["rep"].unique()[0])
                curr_sample = haplo["sample"].unique()[0]
                
                conditions = (true_tracts["ind"].astype(str) == curr_sample) & (true_tracts["rep"].astype(int) == curr_rep)
                curr_truth_tracts = true_tracts[conditions]
                true_tracts.drop(true_tracts[conditions].index)
            
                truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
            
                inferred_tracts = haplo[["chrom", "tractstart", "tractend"]]
                inferred_tracts.columns = ["chrom", "start", "end"]
                truth_tracts.columns = ["chrom", "start", "end"]

                truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
                inferred_tracts = pybedtools.BedTool.from_dataframe(inferred_tracts).sort().merge()

                total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
                total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
                true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

                all_total_inferred_tracts = all_total_inferred_tracts + total_inferred_tracts
                all_total_true_tracts = all_total_true_tracts + total_true_tracts
                all_true_positives = all_true_positives + true_positives
                        
    
    for i, tract in true_tracts.iterrows():
        curr_truth_tracts = tract.to_frame().T
        truth_tracts = curr_truth_tracts[["chr", "start", "end"]]
        truth_tracts.columns = ["chrom", "start", "end"]
        truth_tracts = pybedtools.BedTool.from_dataframe(truth_tracts).sort().merge()
        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        all_total_true_tracts = all_total_true_tracts + total_true_tracts


    if float(all_total_inferred_tracts) == 0: precision = np.nan
    else: precision = all_true_positives / float(all_total_inferred_tracts) * 100
    if float(all_total_true_tracts) == 0: recall = np.nan
    else: recall = all_true_positives / float(all_total_true_tracts) * 100

    return precision, recall

def split_dfs_sample_replicate(df, save_replicates=False):
    """
    Description:
        Helper function for partition the merged data for all replicates and all samples (
        creating complete DataFrames at first and afterwards split them is, of course, redundant work;
        it would be certain beneficial to change the workflow a bit - especially for working on the cluster etc.)

    Arguments:
        df DataFrame: list of dataframes containing true introgresssed tracts.
   
    Returns:
        new_dfs list: list of DataFrames containing information about one sample of one replicate
   
    """
   
    replicates = [v for k, v in df.groupby('rep')]

    if save_replicates == True:
        for i,replicate in enumerate(replicates):
            replicate.to_csv('replicates' + str(i) + ".csv", sep="\t", header=False, index=False)

    new_dfs = []
    for replicate in replicates:
        all_samples = [v for k, v in replicate.groupby('sample')]
        new_samples = []
        for sample in all_samples:
            all_haplos = [v for k, v in sample.groupby('haplo')]
            new_samples.append(all_haplos)
        new_dfs.append(new_samples)

    return new_dfs
