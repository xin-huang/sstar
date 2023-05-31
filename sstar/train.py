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
import numpy as np
from multiprocessing import Process, Queue
import allel
from sstar import stats
from sstar import preprocess
import os


def train(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed=None, train_archie=False):
    """
    """

    # simulate data
    _simulation_manager(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, seed)

    if train_archie:
        #_train_archie()
        _train_archie(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir)
    else:
        _train_sstar()


def _train_archie(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, scikit=True, statsmodels=True, drop_dynamic_cols=True):	

    outputfolder = output_dir

    #set filenames for individuals, reference and target
    ref_ind_file = "new_sim.ref.ind.list"
    tgt_ind_file  = "new_sim.tgt.ind.list"

    create_ref_tgt_file(nref, ntgt, ref_ind_file, tgt_ind_file)

    anc_allele_file = None

    #set window length and stepsize
    win_len = 50000
    win_step = 50000
   
    #I think these parameters are NOT necessary for ArchIE - just retained for the signature of preprocess.process_data
    match_bonus = 1
    max_mismatch = 1
    mismatch_penalty = 1

    process_archie = True

    #tracts with a proportion between not_archaic and archaic are labeled as ambiguous (in _label)
    archaic_prop = 0.7
    not_archaic_prop = 0.3

    true_tracts = []

    true_tracts_labeled = []

    features = []

    file_names = []

    replicate_counter = 0

    #reading of data, preprocessing - i.e., calculating statistics -, and obtaining & labeling of true tracts
    for replicate, file in enumerate(os.listdir(output_dir)):
        if file.endswith(".vcf"):

            filename = os.path.splitext(file)[0]


            feature_file = os.path.splitext(file)[0]+'.features'
            #computation of statistics
            preprocess.process_data(os.path.join(output_dir,file), ref_ind_file, tgt_ind_file, anc_allele_file, os.path.join(output_dir,feature_file), win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, process_archie)


            true_tract = os.path.splitext(file)[0]+'.true.tracts.bed'
            true_tracts.append ( pd.read_csv(os.path.join(output_dir, true_tract), sep="\t", header=None, names=['chr', 'start', 'end', 'hap', 'ind']) )

            
            #labeling of true tracts
            true_tract_labeled = _label(os.path.join(output_dir, true_tract), archaic_prop, not_archaic_prop, seq_len)
            if true_tract_labeled is not None:
                true_tract_labeled["rep"] = replicate_counter

            true_tracts_labeled.append(true_tract_labeled)

            feature = pd.read_csv(os.path.join(output_dir, feature_file), sep="\t")
            feature["rep"] = replicate_counter

            features.append(feature)
            file_names.append(filename)
            replicate_counter = replicate_counter + 1


    feature_df_labeleds = []

    #Labeling of features
    for i, feature_df in enumerate(features):

        feature_df_labeled = label_feature_df_archie(feature_df, true_tracts_labeled[i])

        feature_df_labeleds.append(feature_df_labeled)
    
    #create one big training dataframe
    train_df = pd.concat(feature_df_labeleds)

    #train_df.to_csv(os.path.join(outputfolder, "features_full.csv"))

    #drop all unneccesary columns
    train_df.drop(['rep', 'chrom', 'start', 'end', 'sample', 'interval', 'overlap', 'overlap_percentage', "label_one_1", "label_one_2", "label_one_3", "haplo"], axis=1, inplace=True, errors='ignore')

    #drop_dynamic_cols indicate whether non-fixed size features should be dropped
    if drop_dynamic_cols == True:
        dynamic_cols = [col for col in train_df.columns if ('-ton' in col or col.startswith("pairwised_dist"))]
        train_df.drop(dynamic_cols, axis=1, inplace = True, errors='ignore')


    train_df.to_csv(os.path.join(outputfolder, "features_final.csv"))

    #start training
    scikit_file = output_prefix + ".scikit.pickle"
    statsmodels_file = output_prefix + ".statsmodels.pickle"

    #call training functions

    train_statsmodels(train_df, statsmodels_file)

    train_scikit(train_df, scikit_file)


def create_ref_tgt_file(nref, ntgt, ref_ind_file, tgt_ind_file, identifier="tsk_"):
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



def train_statsmodels(train_df, save_filename):
    """
    Description:
        Function for training of the statsmodels logistic classification.

    Arguments:
        train_df DataFrame: Training data
        save_filename str: filename for output model
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    import numpy as np

    sm_data_exog = train_df.copy()
    sm_data_exog.drop(["label"], axis=1, inplace=True)
    sm_data_exog.replace(np.nan, 0, inplace=True)
    sm_data_exog = sm.add_constant(sm_data_exog, prepend=False)

    sm_data_endog = train_df["label"]

    glm_binom = sm.GLM(sm_data_endog.astype(int), sm_data_exog.astype(float),family=sm.families.Binomial())
    result = glm_binom.fit()

    result.save(save_filename)

    


def train_scikit(train_df, save_filename):
    """
    Description:
        Function for training of the Scikit logistic classification.

    Arguments:
        train_df DataFrame: Training data
        save_filename str: filename for output model
    """

    from sklearn.linear_model import LogisticRegression

    y = train_df['label']
    X = train_df.drop(["label"], axis=1, inplace=False)

    X.replace(np.nan, 0, inplace=True)
    X.replace(pd.NA, 0, inplace=True)

    model = LogisticRegression(solver="newton-cg", penalty=None, max_iter=1000)

    model.fit(X,y.astype(int))

    import pickle

    pickle.dump(model, open(save_filename, "wb"))


def _train_sstar():
    pass




def label_feature_df_archie(feature_df, true_tract_list, discard_ambiguous=True, replicates = True):
    """
    Description:
        Label data for training

    Arguments:
        feature_df DataFrame: (unlabeled) Training data
        true_tract_list DataFrame: true tracts to use for labeling
        discard_ambiguous Boolean: discard tracts classified as ambiguous 
    """

    feature_df['label'] = 0
   
    if true_tract_list is None:
        return feature_df
    
    #print("true tract list")
    #print(true_tract_list)    

    #print("feature_df")
    #print(feature_df)
    true_tract_list["hap"] = true_tract_list["hap"].str.replace("hap_", "")
 
    for ie, entry in enumerate(true_tract_list.values.tolist()):

        haplo = int(entry[0])
        ind = entry[1]
        tract_label = entry[4]
        if replicates == True:
            replicate = int(entry[5])
        


        if replicates == True:
            conditions = (feature_df["sample"] == ind) & (feature_df["haplo"] == haplo) & (feature_df["rep"] == replicate)
        else:
            conditions = (feature_df["sample"] == ind) & (feature_df["haplo"] == haplo) 

        if tract_label[0] == 1:
            
            feature_df.loc[conditions, "label"] = 1
        
        if tract_label[2] == 1:
            
            if discard_ambiguous == True:
                 feature_df = feature_df.drop(feature_df[conditions].index)
        
    return feature_df




def label_feature_df(feature_df, true_tract_list, only_above_threshold=False, discard_ambiguous=True, replicates=False):
    """
    Description:
        Label data for training and compute overlap - NOT USED

    Arguments:
        feature_df DataFrame: (unlabeled) Training data
        true_tract_list DataFrame: true tracts to use for labeling
        only_above_threshold Boolean: only retain 
        discard_ambiguous Boolean: discard tracts classified as ambiguous 
    """

    import numpy as np
    vectorize_overlap = np.vectorize(getOverlap_features_tracts_np)
    vectorize_percentage_overlap = np.vectorize(getOverlap_features_tracts_percentage_np)

    feature_df['label'] = 0
    feature_df['overlap'] = 0
    feature_df['overlap_percentage'] = 0

    feature_df['label_one_1'] = 0
    feature_df['label_one_2'] = 0
    feature_df['label_one_3'] = 0

    label_index = feature_df.columns.get_loc("label")
    label_one_index = feature_df.columns.get_loc("label_one_1")
    label_two_index = feature_df.columns.get_loc("label_one_2")
    label_three_index = feature_df.columns.get_loc("label_one_3")


    overlap_index = feature_df.columns.get_loc("overlap")
    overlap_percentage_index = feature_df.columns.get_loc("overlap_percentage")
    start_index = feature_df.columns.get_loc("start")
    end_index = feature_df.columns.get_loc("end")
    haplo_index = feature_df.columns.get_loc("haplo")
    sample_index = feature_df.columns.get_loc("sample")
    if replicates == True:
        rep_index = feature_df.columns.get_loc("rep")

    feature_array = feature_df.to_numpy()


    true_tract_list["hap"] = true_tract_list["hap"].str.replace("hap_", "")
    

    for ie, entry in enumerate(true_tract_list.values.tolist()):
        start = entry[1]
        end = entry[2]

        haplo = int(entry[3])
        ind = entry[4]
        if replicates == True:
            replicate = int(entry[5])
        
        if replicates == True:
            conditions=np.where((feature_array[:, sample_index] == ind) & ( feature_array[:, haplo_index]  == haplo) & ( feature_array[:, rep_index]  == replicate))
        else:
            conditions=np.where((feature_array[:, sample_index] == ind) & ( feature_array[:, haplo_index]  == haplo))
        
        feature_array[conditions, overlap_index] = vectorize_overlap(feature_array[conditions, overlap_index], start, end, feature_array[conditions,start_index], feature_array[conditions,end_index]) 
        
        feature_array[conditions, overlap_percentage_index] = vectorize_percentage_overlap(feature_array[conditions, overlap_percentage_index], start, end, feature_array[conditions,start_index], feature_array[conditions,end_index]) 

        feature_array[np.where((feature_array[:, overlap_index] > 0)), label_index] = 1
        

    feature_array[np.where((feature_array[:, overlap_percentage_index] >= 0.7)), label_one_index] = 1
    feature_array[np.where((feature_array[:, overlap_percentage_index] < 0.7) & (feature_array[:, overlap_percentage_index] >= 0.3)), label_two_index] = 1
    feature_array[np.where((feature_array[:, overlap_percentage_index]  < 0.3)), label_three_index] = 1

    #back to pandas df for convenience

    feature_df = pd.DataFrame(feature_array, columns = feature_df.columns)

    if only_above_threshold == True:
        feature_df['label'] = feature_df['label_one_1']

    if discard_ambiguous == True:
        feature_df = feature_df.drop(feature_df[feature_df.label_one_2 == 1].index)

    return feature_df

#functions for calculating overlap

def getOverlap_features_tracts_np(startvalue, tract_start, tract_end, feature_start, feature_end):
    return startvalue + max(0, min(tract_end, feature_end) - max(tract_start, feature_start))

def getOverlap_features_tracts_percentage_np(startvalue, tract_start, tract_end, feature_start, feature_end):
    overlap = max(0, min(tract_end, feature_end) - max(tract_start, feature_start))
    percentage = overlap / (feature_end - feature_start)
    percentage = percentage + startvalue
    return percentage


def getOverlap(tract_start, tract_end, feature_start, feature_end):
    return max(0, min(tract_end, feature_end) - max(tract_start, feature_start))

def getOverlapPercentage(tract_interval, feature_interval):
    overlap = max(0, min(tract_interval[1], feature_interval[1]) - max(tract_interval[0], feature_interval[0]))
    percentage = overlap / (feature_interval[1] - feature_interval[0])
    return percentage

def getOverlap_pd(tract_int, feature_int):
    feature_start = feature_int.left
    feature_end = feature_int.right
    tract_start = tract_int.left
    tract_end = tract_int.right
    return max(0, min(feature_end, tract_end) - max(feature_start, tract_start))

def getOverlapPercentage_pd(tract_int, feature_int):
    feature_start = feature_int.left
    feature_end = feature_int.right
    tract_start = tract_int.left
    tract_end = tract_int.right
    overlap = max(0, min(feature_end, tract_end) - max(feature_start, tract_start))
    percentage = overlap / (feature_end - feature_start)
    return percentage


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

    global df

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
        for n in sorted(true_tracts.keys()):
            true_tracts[n].sort(key=lambda x:(x[0], x[1], x[2]))
            df2 = pd.DataFrame(true_tracts[n], columns=['chr', 'start', 'end', 'hap', 'ind'])
            df = pd.concat([df, df2])

        df.drop_duplicates(keep='first').to_csv(output_dir+'/'+output_prefix+f'{rep}.true.tracts.bed', sep="\t", header=False, index=False)



        out_queue.put(rep)


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
            if i.left > t.interval[0]: continue
            if i.right <= t.interval[0]: break
            for n in tracts.keys():
                if t.is_descendant(n, i.node): tracts[n].append([1, int(i.left), int(i.right), f'hap_{int(n%2)}', f'tsk_{ts.node(n).individual}'])

    return tracts


def _label(tracts, archaic_prop, not_archaic_prop, seq_len):
    """
    Description:
        Helper function to label a fragment as 'introgressed', 'not introgressed', or 'ambiguous'.

    Arguments:
        tracts str: Name of the file containing ground truth introgressed fragments.
        archaic_prop float: Threshold to label a fragment as 'introgressed'.
        not_archaic_prop float: Threshold to label a fragment as 'not introgressed'.
        seq_len int: Length of the fragment.
    """
    def _add_label(row, archaic_prop, not_archaic_prop):
        if row['prop'] > archaic_prop: return [1,0,0]
        elif row['prop'] < not_archaic_prop: return [0,1,0]
        else: return [0,0,1]

    try:
        df = pd.read_csv(tracts, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return None
    
    df = pd.read_csv(tracts, sep="\t", header=None)

    df.columns = ['chr', 'start', 'end', 'hap', 'ind']
    df['len'] = df['end'] - df['start']
    df = df.groupby(by=['hap', 'ind'])['len'].sum().reset_index()
    df['prop'] = df['len'] / seq_len
    df['label'] = df.apply(lambda row: _add_label(row, archaic_prop, not_archaic_prop), axis=1)

    #sep="\t", header=None, names=['chr', 'start', 'end', 'hap', 'ind']
    return df




def train_parameters_archienew_fin(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir):

    preprocess.store_global_parameters(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir)

    train(demo_model_file, nrep, nref, ntgt, ref_id, tgt_id, src_id, seq_len, mut_rate, rec_rate, thread, output_prefix, output_dir, train_archie=True)


if __name__ == '__main__':
    outputfolder = "example_output"
    train("./examples/models/BonoboGhost_4K19.yaml", 100, 50, 50, 'Western', 'Bonobo', 'Ghost', 50000, 1e-8, 1e-8, 2, 'testa', outputfolder, train_archie=True)

