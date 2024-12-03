# identification-of-TCR-sequences-in-GCA

## Author Contribution Statement
All code in this repository was sourced from [KNN_classification_of_T_cell_receptor_sequences_associated_with_giant_cell_arteritis](https://github.com/annaweber209/KNN_classification_of_T_cell_receptor_sequences_associated_with_giant_cell_arteritis), originally written by [annaweber209](https://github.com/annaweber209), and modified by [mlzski](https://github.com/mlzski) for this article: **Identification of clonally expanded T-cell receptor sequences in giant cell arteritis - REF HERE**.

## Installation

The library itself has few dependencies (see [setup.py](setup.py)) with loose requirements. 

Create a virtual environment and install dependencies

```console
python -m venv --system-site-packages venv
source venv/bin/activate
pip install -r requirements.txt
```
Install in editable mode for development:
```console
pip install -e .
```


## Data structure
If you want to apply the code to your own data, create two .csv files of TCR sequences containing the following columns:
```console
cloneCount, cloneFraction, cdr3_b_aa, v_b_gene, j_b_gene, Name, Label
```
where cdr3_b_aa is the amino acid sequence of the CDR3 region, v_b_gene and j_b_gene are the names of the V and J segments in IMGT format (e.g. TRBV5-6\*01), Name is the Name/ID of the patient that this TCR was found in and Label is the label of that patient (i.e. 0 = control, 1 = disease).



# Example usage
## Train an ensemble of K-NN models
### Get pairwise distances
All the parameters and hyperparameters can be specified in a single parameter file, like the example shown in parameter_files.

```console
python3 get_sparse_dist_matrix.py \
parameter_files/dist_matrix.json
```

The output is a file containing all pairwise distances of the input TCRs. Note that the computation is memory-intensive, but can be broken into smaller batches. Default batch size is set to 5000, but if this fills up your memory you can lower it for longer but safe computation. In general, it is advisable to run this code either on small datasets or on a virtual machine/cluster. 

### Classify TCRs
To predict the probability of each TCR to belong to the diseased class, run the below code. 
Again, all parameters and hyperparameters of the model can be specified in a single parameter file. Here you can set how many independent experiments to run, how many positive and negative TCRs to draw from your input data (recommended to keep this <10,000), which number of neighbours k to use for the K-NN, how many patient-wise cross-validation folds to use and whether to randomize labels for a bootstrapping experiment.

```console
python3 clustering_from_saved_dist.py \
parameter_files/classifier.json
```

This code outputs a lot of files (number of rounds times number of cv folds) with columns

```console
v_b_gene,j_b_gene,cdr3_b_aa,GCA,Name,labels,predictions
```

### Cluster TCRs
Finally, we provide code that identifies the k nearest neighbors of a set of TCRs (query) within the full TCR input (all controls and disease TCRs). Note that it is not advised to query more than 10,000 TCRs at once.

```console
python3 get_NN.py \
parameter_files/get_NN.json
```


## Data Handling
In utils you find some functions for convenient handling of TCR sequences. 
To generate full sequences of TCRs from CDR3 sequence and V and J segment names, the `cdr3_to_full_seq.py` script can be used. The script relies on the user having downloaded a fasta files containing the Names of V and J segments with their respecive sequences called `V_segment_sequences.fasta` and `J_segment_sequences.fasta`. These can be downloaded from IMGT.org. With the `align.py` file CDR3 sequences can be aligned according to the IMGT scheme, e.g. for generation of sequence logos. 








