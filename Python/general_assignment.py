import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from matplotlib import cm
import argparse

parser = argparse.ArgumentParser(description='Parser for Summary Results')

parser.add_argument('--htodemul_assignment',  help='CSV file with Hashtag assignment results from HTOdemul', default="-")
parser.add_argument('--multiseq_results',  help='CSV file with classification results from MULTI-seq',default="-")
parser.add_argument('--hashed_drop_results',  help='CSV file with classification results from Hashed Drop',default="-")
parser.add_argument('--demuxem_results',  help='CSV file with classification results from DemuxEM',default="-")
parser.add_argument('--hashsolo_results',  help='CSV file with classification results from Hash Solo',default="-")

parser.add_argument('--output_assignment',  help='name for the output file',default="general_assignment.csv")

args = parser.parse_args()



total_classification = pd.DataFrame()

if(args.hashsolo_results != "-"): 
    hashsolo_results = pd.read_csv(args.hashsolo_results)
    #Processing Hash Solo results
    hashsolo_hash = hashsolo_results.copy()
    del hashsolo_hash["most_likely_hypothesis"]
    del hashsolo_hash["cluster_feature"]
    del hashsolo_hash["negative_hypothesis_probability"]
    del hashsolo_hash["singlet_hypothesis_probability"]
    del hashsolo_hash["doublet_hypothesis_probability"]
    hashsolo_hash.columns = ["Barcode","Assignment-hashsolo"]
    total_classification.append(hashsolo_hash)

if(args.demuxem_results != "-"):  
    demuxem_results = pd.read_csv(args.demuxem_results)
    #Processing DemuxEM results
    demuxem_hashes = demuxem_results.copy()
    del demuxem_hashes['n_genes']
    del demuxem_hashes['n_counts']
    del demuxem_hashes['passed_qc']
    del demuxem_hashes['demux_type']
    demuxem_hashes.columns = ["Barcode_demux","Assignment-demuxem"]
    demuxem_hashes[['Assignment-demuxem']] = demuxem_hashes[['Assignment-demuxem']].fillna('Negative')
    demuxem_hashes.loc[demuxem_hashes['Assignment-demuxem'].str.contains(','), 'Assignment-demuxem'] = 'Doublet'
    if(args.hashsolo_results != "-"):
        #only delete barcode section if hash solo is out of the picture
        del demuxem_hashes["Barcode_demux"]
        total_classification["Assignment-demuxem"]=demuxem_hashes
    else:
        total_classification = pd.concat([total_classification,demuxem_hashes],axis=1)
        


if(args.htodemul_assignment != "-"):
    htodemul_hashes = pd.read_csv(args.htodemul_assignment) 
    #Processing HTODemux results
    htodemul_hashes.columns = ["Barcode_htoDemux","Assignment-htoDemux"]
    htodemul_hashes.loc[htodemul_hashes['Assignment-htoDemux'].str.contains('_'), 'Assignment-htoDemux'] = 'Doublet'
    if(args.hashsolo_results == "-" and args.demuxem_results == "-"):
        total_classification = pd.concat([total_classification,htodemul_hashes],axis=1)
    else:
        del htodemul_hashes["Barcode_htoDemux"]
        total_classification["Assignment-htoDemux"] = htodemul_hashes



if(args.multiseq_results != "-"):
    multiseq_results = pd.read_csv(args.multiseq_results)
    #Processing MULTI-seq results
    multiseq_hashes = multiseq_results.copy()
    multiseq_hashes.columns = ["Barcode_multi","Assignment-multiseq"]
    if(args.hashsolo_results == "-" and args.demuxem_results == "-" and args.multiseq_results  == "-"):
        total_classification = pd.concat([total_classification,multiseq_hashes],axis=1)
    else:
        del multiseq_hashes["Barcode_multi"]
        total_classification["Assignment-multiseq"] = multiseq_hashes

if(args.hashed_drop_results != "-"):   
    hasheddrops_results = pd.read_csv(args.hashed_drop_results)
    #Processing Hashed Drop Results
    #Tranform singlet - doublet column into 1s and 0s
    hasheddrops_results["Doublet"] = hasheddrops_results["Doublet"].astype(int)
    hasheddrops_results["Confident"] = hasheddrops_results["Confident"].astype(int)
    #List of conditions to be applied between singlet and doublet columns, so that we can determine which one is which
    condlist = [(hasheddrops_results['Doublet']>hasheddrops_results['Confident']),(hasheddrops_results['Doublet']<hasheddrops_results['Confident']),(hasheddrops_results['Doublet']==hasheddrops_results['Confident'])]
    #List with the classification labels
    choicelist = ["Doublet", "Singlet","Negative"]
    hasheddrops_results['HashedDrops'] = np.select(condlist, choicelist, default=0)
    hasheddrops_hashes = hasheddrops_results.copy()
    del hasheddrops_hashes['Total']
    del hasheddrops_hashes['Second']
    del hasheddrops_hashes['LogFC']
    del hasheddrops_hashes['LogFC2']
    del hasheddrops_hashes['Doublet']
    del hasheddrops_hashes['Confident']
    del hasheddrops_hashes['HashedDrops']
    hasheddrops_hashes.columns = ["Barcode_hashed","Assignment-hashed"]
    if(args.hashsolo_results == "-" and args.demuxem_results == "-" and args.multiseq_results  == "-" and args.multiseq_results== "-"):
        total_classification = pd.concat([total_classification,hasheddrops_hashes],axis=1)
    else:
        del hasheddrops_hashes["Barcode_hashed"]
        total_classification["Assignment-hashed"] = hasheddrops_hashes



total_classification.to_csv("general_assignment.csv")