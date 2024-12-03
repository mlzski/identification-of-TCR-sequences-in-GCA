import numpy as np
import pandas as pd
import argparse

#MiXCR does not always aggregate all sequences with identical CDR3 sequence
#or even identical CDR3 sequence and V segment into the same clone
#this script identifies these instances and merges clones
#by default sequences with identical CDR3 sequence are merged
#the --V option lets you merge sequences with identical CDR3 and V segment

parser = argparse.ArgumentParser()
parser.add_argument(
    'data_filepath', type=str,
    help='Path to exported clones data file.'
)
parser.add_argument(
    'output_filepath', type=str,
    help='Path to folder where output should be stored.'
)
parser.add_argument(
    '--V', action='store_true',
    help='Use this option to cluster by CDR3 AND V segment identity.'
)


def main(data_filepath, output_filepath, V):
    # Read in clones data from MiXCR output file
    data = pd.read_csv(data_filepath, delimiter='\t')
    if V:
        data['identifier'] = data['aaSeqCDR3']+data['bestVHit']
    else:
        data['identifier'] = data['nSeqCDR3']
    data_new = pd.DataFrame(columns=data.columns)

    for i in np.arange(len(data)):
        cdr3 = data.iloc[i]['nSeqCDR3']
        vseg = data.iloc[i]['bestVHit']
        if V:
            condition = cdr3+vseg in data_new['identifier'].tolist()
            id = cdr3+vseg
        else:
            condition = cdr3 in data_new['nSeqCDR3'].tolist()
            id = cdr3
        
        if condition:
            # add x to cloneCount of that row in data_new
            clonecount = int(data_new.loc[data_new['identifier'] == id]['cloneCount'])
            new_clonecount = clonecount + int(data.iloc[i]['cloneCount'])
            data_new['cloneCount'] = np.where(
                data_new['identifier'] == id,
                new_clonecount,
                data_new['cloneCount']
                )
        else:
            data_new = data_new.append(data.iloc[i])

    print('Read data from', data_filepath)
    print('Total Number of Counts', sum(data['cloneCount']))
    print('Total Number of Individual Clones before', len(data))
    print('Total Number of Individual Clones after', len(data_new))
    assert sum(data['cloneCount']) == sum(data_new['cloneCount'])
    data_new = data_new.drop(columns=['identifier'])
    data_new.to_csv(output_filepath, index=None, sep='\t', mode='w')


if __name__ == '__main__':
    # parse arguments
    args = parser.parse_args()
    main(args.data_filepath, args.output_filepath, args.V)