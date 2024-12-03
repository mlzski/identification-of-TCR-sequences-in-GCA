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
import sys
from scipy.sparse import vstack, save_npz, load_npz
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

    X = tr.rw_beta # sparse data distances matrix
    return X


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
    all = process(all_controls_cleaned, all_GCA_cleaned)
    nr_batches = int(len(all)/params['batchsize'])

    for i, chunk in enumerate(np.array_split(all, nr_batches)):
        logger.info(f"batch {i} of {nr_batches}")
        X = get_distances(all,chunk, params)
        if i == 0:
            distances = X
        else:
            distances = vstack([distances, X])
    save_npz(os.path.join(params['outfolder'], params['outfile']+"_distance_matrix.npz"), distances)
    all.to_csv(os.path.join(params['outfolder'], params['outfile']+"raw_data.csv"), index=False)

if __name__ == '__main__':
     # parse arguments
     args = parser.parse_args()
     # run the training
     main(
         args.parameter_file
     )