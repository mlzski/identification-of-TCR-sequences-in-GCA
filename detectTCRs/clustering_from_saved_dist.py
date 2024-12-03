import pandas as pd
import os
import numpy as np
from tcrdist.repertoire import TCRrep
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from scipy.sparse import vstack, save_npz, load_npz
from sklearn import preprocessing
import logging 
import json 

import sys
import argparse 


def get_subsample(data, params):
    data = data.drop_duplicates(subset='cdr3_b_aa')
    data = data[["v_b_gene", "j_b_gene", "cdr3_b_aa", "GCA", "Name"]]
    data['cdr3_b_aa'] = data['cdr3_b_aa'].apply(lambda x: x.upper())

    all_GCA = data[data['GCA'] == 1]
    all_controls =  data[data['GCA'] == 0]
    positives_sample = all_GCA.sample(n=params['nr_positives']) #random subset
    negatives_sample = all_controls.sample(n=len(positives_sample))

    output_data = pd.concat([negatives_sample, positives_sample])
    
    return output_data

def dist_unsparse(dist_matrix, subsample, params):
    # select sampled rows from dist matrix
    indices = list(subsample.index)
    X = dist_matrix[indices, :]
    X = X[:, list(subsample.index)]
    X = np.absolute(X.toarray()) # dense training data distances matrix
    Y = subsample['GCA'] # training data labels
    # set all zeros to 100 (max distance) and normalize
    X[X == 0] = params['max_distance']
    X = pd.DataFrame(preprocessing.normalize(X)) 
    return X, Y
        
def get_split(X,Y,params):
    # 3) Split into CV; ensure all TCRs from one subject are in the same group
    cv_folds=params['cv_folds']
    kfolds = StratifiedGroupKFold(n_splits=cv_folds, shuffle=True, random_state=42)
    kfolds.get_n_splits(X, Y)
    return kfolds

def fit_model(X_train,Y_train,X_test,Y_test,params):
    # 4) fit K-NN to train distances, prediction for val distances
    model_fold = KNeighborsClassifier(n_neighbors=params['neighbors'], metric='precomputed', weights='distance')
    model_fold.fit(X_train, Y_train.values.ravel())
    predictions = model_fold.predict_proba(X_test)[:,1] #take only probability for predicting 1
    assert len(Y_test) == len(predictions)
    return predictions

parser = argparse.ArgumentParser()
parser.add_argument(
    'parameter_file', type=str,
    help='Path to the parameter file.'
)

def main(
    parameter_file
):

    # setup logging
    logging.basicConfig(stream=sys.stdout)
    logger = logging.getLogger(f'{parameter_file}')
    logger.setLevel(logging.DEBUG)
    with open(parameter_file, "r") as jsonfile:
         params = json.load(jsonfile)

    # Load data and distance matrix
    all = pd.read_csv(os.path.join(params['infolder'], params['data_file']+"raw_data.csv"))
    loaded_dist = load_npz(os.path.join(params['infolder'], params['distance_file']+"_distance_matrix.npz"))

    for round in range(params['rounds']):
        run_id = params['id_start']+round
        subsample = get_subsample(all, params)
        logger.info(f"Subsample of lenght {len(subsample)}")
        logger.info(f"{len(subsample[subsample['GCA'] == 1])} positive data")
        logger.info(f"{len(subsample[subsample['GCA'] == 0])} negative data")

        X,Y = dist_unsparse(loaded_dist, subsample, params)
        kfolds = get_split(X,Y,params)

        fold = 0
        for train_ind, test_ind in kfolds.split(X, Y, subsample['Name']):
            tcrs_test = subsample.iloc[test_ind]
            X_train, X_test = X.iloc[train_ind, train_ind], X.iloc[test_ind, train_ind]
            Y_train, Y_test = Y.iloc[train_ind], Y.iloc[test_ind]
            if params['randomize'] == 'True':
                Y_train = Y_train.sample(frac=1)
            predictions = fit_model(X_train,Y_train,X_test,Y_test,params)

            # Write out relevant info
            tcrs_test['labels'] = Y_test.values.ravel()
            tcrs_test['predictions'] = predictions
            tcrs_test.to_csv(os.path.join(params['outfolder'], str(run_id)+'_'+str(fold)+'_k'+str(params['neighbors'])+'_predictions.csv'), index=False)

            fold+=1
        logger.info(f"Finished with round {round}")

if __name__ == '__main__':
     # parse arguments
     args = parser.parse_args()
     # run the training
     main(
         args.parameter_file
     )