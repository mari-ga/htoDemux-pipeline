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

parser.add_argument('--output_assignment',  help='name for the output file',default="general_classification.csv")

args = parser.parse_args()

htodemux = False
multiseq = False
demuxem = False
hasheddrops = False
hashsolo = False




total_classification = pd.DataFrame()

#if(args.hashsolo_results != "-"): 
if(args.hashsolo_results != "-" ):
    #Flag so that I can know which tools produced results and which ones didn;t
    hashsolo = True
    hashsolo_results = pd.read_csv(args.hashsolo_results)
    #Processing Hash Solo results
    del hashsolo_results["cluster_feature"]
    del hashsolo_results["negative_hypothesis_probability"]
    del hashsolo_results["singlet_hypothesis_probability"]
    del hashsolo_results["doublet_hypothesis_probability"]
    del hashsolo_results["Classification"]
    hashsolo_results.columns = ["Barcode-total","Classification-HashSolo"]
    #List of conditions to be applied between singlet and doublet columns, so that we can determine which one is which
    condlist_solo = [(hashsolo_results['Classification-HashSolo']==1),
                 (hashsolo_results['Classification-HashSolo']==0),
                 (hashsolo_results['Classification-HashSolo']==2)]
    #List with the classification labels
    choicelist_solo = ["Singlet", "Negative","Doublet"]
    hashsolo_results['Classification-HashSolo'] = np.select(condlist_solo, choicelist_solo)             
    total_classification= pd.concat([total_classification,hashsolo_results],axis=1)
    

#if(args.demuxem_results != "-"):  
if(args.demuxem_results != "-"):
    demuxem = True
    demuxem_results = pd.read_csv(args.demuxem_results)
    #Processing DemuxEM results
    del demuxem_results['n_genes']
    del demuxem_results['n_counts']
    del demuxem_results['passed_qc']
    del demuxem_results['assignment']
    demuxem_results.columns = ["Barcode_demux","Classification-DemuxEM"]
    demuxem_results.loc[demuxem_results['Classification-DemuxEM'].str.contains('doublet'), 'Classification-DemuxEM'] = 'Doublet'
    demuxem_results.loc[demuxem_results['Classification-DemuxEM'].str.contains('singlet'), 'Classification-DemuxEM'] = 'Singlet'
    demuxem_results.loc[demuxem_results['Classification-DemuxEM'].str.contains('unknown'), 'Classification-DemuxEM'] = 'Negative'

    if(hashsolo == True):
        #only delete barcode section if hash solo is out of the picture
        del demuxem_results["Barcode_demux"]
        total_classification["Classification-DemuxEM"]=demuxem_results
    else:
        total_classification = pd.concat([total_classification,demuxem_results],axis=1)
        


if(args.htodemul_classification != "-"):
    htodemux = True
    htodemux_classification = pd.read_csv(args.htodemul_classification) 
    #Processing HTODemux results
    htodemux_classification.columns = ["Barcode_htoDemux","Classification-HTODemux"]
    if(hashsolo == False and demuxem == False):
        total_classification = pd.concat([total_classification,htodemux_classification],axis=1)
    else:
        del htodemux_classification["Barcode_htoDemux"]
        total_classification["Classification-HTODemux"] = htodemux_classification



if(args.multiseq_results != "-"):
    multiseq = True
    multiseq_results = pd.read_csv(args.multiseq_results)
    #Processing MULTI-seq results
    multiseq_results.columns = ["Barcode_multi","Classification-MultiSeq"]
    multiseq_results.loc[multiseq_results['Classification-MultiSeq'].str.contains('Hash'), 'Classification-MultiSeq'] = 'Singlet'
    if(hashsolo == False  and demuxem == False and htodemux  == False):
        total_classification = pd.concat([total_classification,multiseq_results],axis=1)
    else:
        del multiseq_results["Barcode_multi"]
        total_classification["Classification-MultiSeq"] = multiseq_results

#if(args.hashed_drop_results != "-"):   
if(args.hashed_drop_results != "-"):
    hasheddrops = True
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
    if(hashsolo == False  and demuxem == False and htodemux  == False and multiseq == False):
        total_classification = pd.concat([total_classification,hasheddrops_results],axis=1)
    else:
        del hasheddrops_results["Barcode_hashed"]
        total_classification["Classification-HashedDrops"] = hasheddrops_results




total_classification.to_csv(args.output_assignment)