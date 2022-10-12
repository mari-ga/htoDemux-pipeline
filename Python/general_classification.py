import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from matplotlib import cm
import argparse

parser = argparse.ArgumentParser(description='Parser for Summary Results')

parser.add_argument('--htodemul_classification',  help='CSV file with Hashtag classification results from HTOdemul', default="-")
parser.add_argument('--multiseq_results',  help='CSV file with classification results from MULTI-seq',default="-")
parser.add_argument('--hashed_drop_results',  help='CSV file with classification results from Hashed Drop',default="-")
parser.add_argument('--demuxem_results',  help='CSV file with classification results from DemuxEM',default="-")
parser.add_argument('--hashsolo_results',  help='CSV file with classification results from Hash Solo',default="-")
parser.add_argument('--solo_results',  help='CSV file with prediction results from Solo',default="-")

parser.add_argument('--output_assignment',  help='name for the output file',default="general_classification.csv")

args = parser.parse_args()



total_classification = pd.DataFrame()

if(args.hashsolo_results != "-"): 
    hashsolo_results = pd.read_csv(args.hashsolo_results)
    #Processing Hash Solo results
    hashsolo_hash = hashsolo_results.copy()
    del hashsolo_hash["cluster_feature"]
    del hashsolo_hash["negative_hypothesis_probability"]
    del hashsolo_hash["singlet_hypothesis_probability"]
    del hashsolo_hash["doublet_hypothesis_probability"]
    del hashsolo_results["Classification"]
    hashsolo_hash.columns = ["Barcode-total","Classification-HashSolo"]
    #List of conditions to be applied between singlet and doublet columns, so that we can determine which one is which
    condlist_solo = [(hashsolo_results['Classification-HashSolo']==1),
                 (hashsolo_results['Classification-HashSolo']==0),
                 (hashsolo_results['Classification-HashSolo']==2)]
    #List with the classification labels
    choicelist_solo = ["Singlet", "Negative","Doublet"]
    hashsolo_results['Classification-HashSolo'] = np.select(condlist_solo, choicelist_solo)             
    total_classification= pd.concat([total_classification,hashsolo_hash],axis=1)
    

if(args.demuxem_results != "-"):  
    demuxem_results = pd.read_csv(args.demuxem_results)
    #Processing DemuxEM results
    del demuxem_results['n_genes']
    del demuxem_results['n_counts']
    del demuxem_results['passed_qc']
    del demuxem_results['assignment']
    demuxem_results.columns = ["Barcode_demux","Classification-DemuxEM"]
    demuxem_results.loc[demuxem_results['DemuxEM'].str.contains('doublet'), 'DemuxEM'] = 'Doublet'
    demuxem_results.loc[demuxem_results['DemuxEM'].str.contains('singlet'), 'DemuxEM'] = 'Singlet'
    demuxem_results.loc[demuxem_results['DemuxEM'].str.contains('unknown'), 'DemuxEM'] = 'Negative'

    if(args.hashsolo_results != "-"):
        #only delete barcode section if hash solo is out of the picture
        del demuxem_results["Barcode_demux"]
        total_classification["Classification-DemuxEM"]=demuxem_results
    else:
        total_classification = pd.concat([total_classification,demuxem_results],axis=1)
        


if(args.htodemul_classification != "-"):
    htodemux_classification = pd.read_csv(args.htodemul_classification) 
    #Processing HTODemux results
    htodemux_classification.columns = ["Barcode_htoDemux","Classification-HTODemux"]
    if(args.hashsolo_results == "-" and args.demuxem_results == "-"):
        total_classification = pd.concat([total_classification,htodemux_classification],axis=1)
    else:
        del htodemux_classification["Barcode_htoDemux"]
        total_classification["Classification-HTODemux"] = htodemux_classification



if(args.multiseq_results != "-"):
    multiseq_results = pd.read_csv(args.multiseq_results)
    #Processing MULTI-seq results
    multiseq_results.columns = ["Barcode_multi","Classification-MultiSeq"]
    multiseq_results.loc[multiseq_results['Classification-MultiSeq'].str.contains('Hash'), 'Classification-MultiSeq'] = 'Singlet'
    if(args.hashsolo_results == "-" and args.demuxem_results == "-" and args.htodemul_classification  == "-"):
        total_classification = pd.concat([total_classification,multiseq_results],axis=1)
    else:
        del multiseq_results["Barcode_multi"]
        total_classification["Classification-MultiSeq"] = multiseq_results

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
    del hasheddrops_results['Total']
    del hasheddrops_results['Best']
    del hasheddrops_results['Second']
    del hasheddrops_results['LogFC']
    del hasheddrops_results['LogFC2']
    del hasheddrops_results['Doublet']
    del hasheddrops_results['Confident']
    hasheddrops_results.columns = ["Barcode_hashed","Classification-HashedDrops",]
    if(args.hashsolo_results == "-" and args.demuxem_results == "-" and args.multiseq_results  == "-" and args.multiseq_results== "-"):
        total_classification = pd.concat([total_classification,hasheddrops_results],axis=1)
    else:
        del hasheddrops_results["Barcode_hashed"]
        total_classification["Classification-HashedDrops"] = hasheddrops_results


if(args.solo_results != "-"): 
    solo_results = pd.read_csv(args.solo_results)
    if(args.hashsolo_results == "-" and args.demuxem_results == "-" and args.multiseq_results  == "-" and args.multiseq_results== "-" and args.hasheddrops_hashes == "-"):
        if solo_results.shape[1] == 2:
            solo_results.columns = ["Barcode", "Prediction"]
            total_classification = pd.concat([total_classification,solo_results],axis=1)
        else:
            solo_results.columns = ["Barcode","Doublet", "Singlet"]
            total_classification = pd.concat([total_classification,solo_results],axis=1)
    else:
        if solo_results.shape[1] == 2:
            del solo_results["Barcode"]
            total_classification["Assignment-hashed"] = hasheddrops_hashes
        else:
            total_classification["Doublets"] = hasheddrops_hashes["Doublet"]
            total_classification["Singlets"] = hasheddrops_hashes["Singlet"]
        
   


total_classification.to_csv(args.output_assignment)