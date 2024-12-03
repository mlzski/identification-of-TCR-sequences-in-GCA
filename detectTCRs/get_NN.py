import pandas as pd
import os
import numpy as np
from tcrdist.repertoire import TCRrep
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import (
    auc, average_precision_score, precision_recall_curve, roc_curve
)
from sklearn import preprocessing
import logging 
import json 
import scipy.sparse
import sys
import argparse 


def process(negatives, positives):
    all = pd.concat([negatives, positives], ignore_index=True)
    all = all.drop_duplicates(subset='cdr3_b_aa')
    all = all[["v_b_gene", "j_b_gene", "cdr3_b_aa", "GCA", "Name"]]
    all['cdr3_b_aa'] = all['cdr3_b_aa'].apply(lambda x: x.upper())
    return all

def get_distances(all, to_analyze, params):
    # 2) Calculate distances between all samples (sparse matrix)
    tcrrep_all = TCRrep(cell_df=all,
                        organism='human', chains=['beta'],
                        compute_distances=False,
                        deduplicate=False)
    tcrrep_all.cpus = params['nr_cpus']
    #tcrrep_train.compute_sparse_rect_distances(radius = params['max_distance'], chunk_size = 100)
    df_search = tcrrep_all.clone_df.copy()                 #(3)

    tr = TCRrep(cell_df = to_analyze,                           #(4)
                organism='human', chains=['beta'],
                        compute_distances=False,
                        deduplicate=False)
    tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = df_search, radius = params['max_distance'], chunk_size = 100) #(5)

    # Parse to dense matrix without negative values
    X = tr.rw_beta # dense training data distances matrix
    return X

def output(all, NN_array, params):
    NN_array = pd.DataFrame(NN_array, columns=['TCR1', 'TCR2', 'Distance'])
    NN_array.to_csv(os.path.join(params['outfolder'], 'NN'+str(params['neighbors'])+'.csv'), index=False)
    all.to_csv(os.path.join(params['outfolder'], 'NN_tcrs.csv'), index=False)
    return 

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

    all_controls_cleaned = pd.read_csv(os.path.join(params['infolder'], params['controls']))
    all_GCA_cleaned = pd.read_csv(os.path.join(params['infolder'], params['GCA']))
    tcrs_to_analyze = pd.read_csv(os.path.join(params['infolder'], params['to_analyze']))
    NN_array = []
    
    all = process(all_controls_cleaned, all_GCA_cleaned)
    X = get_distances(all,tcrs_to_analyze, params)
    logger.info(f'Calculated Sparse Distance Matrix')
    scipy.sparse.save_npz('/tmp/sparse_matrix.npz', X)

    for i in range(X.shape[0]):
        row = np.absolute(X[i].toarray())[0]
        row[row == 0] = params['max_distance']
        res = sorted(range(len(row)), key = lambda sub: row[sub])[:params['neighbors']+1]
        for id in res:
            NN_array.append([i, id, row[id]])
    logger.info(f'Finished saving row {i}')

    output(all, NN_array, params)

if __name__ == '__main__':
     # parse arguments
     args = parser.parse_args()
     # run the training
     main(
         args.parameter_file
     )