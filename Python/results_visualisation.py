from ast import arg
import string
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser(description='Parser for Demultiplexing results visualisation ')
parser.add_argument('--htoDemuxResultPath',  help='Path for the csv file containing the results for HTODemux')
parser.add_argument('--multiSeqResultPath',  help='Path for the csv file containing the results for MULTI-seq')
parser.add_argument('--hashedDropsResultPath',  help='Path for the csv file containing the results for Hashed Drops')
parser.add_argument('--graphPath',  help='Path to save graphs')
parser.add_argument('--numberAlgoritms',  help='Number of algorithms that were executed and will be included in the graphs')
args = parser.parse_args()

def pre_processing(path,algorithm):
    results = pd.read_csv(path)
    results.columns = ["hashtag","classification"]
    counts = results['classification'].value_counts()
    results_dict = counts.to_dict()
    results_dict['Algorithm'] = algorithm
    if results_dict.get("Doublet") is None:
        results_dict['Doublet'] = 0
    if results_dict.get("Negative") is None:   
        results_dict['Negative'] = 0

    return results_dict

#The next 2 functions receive dicts as inputs with the results preprocessed so that each one of them produces
#as final result a dict with Singlets, Doublets and Negatives

def multi_seq_visual(results):
    multiseq_dict = dict()
    multi_copy = results.copy()
    del multi_copy['Doublet']
    del multi_copy['Negative']
    del multi_copy['Algorithm']
    #sum add up all singlets 
    multi_singlets = sum(multi_copy.values())
    #multi_singlets = sum((value for key, value in results.items() if key != 'Doublet' or key != 'Negative' or key != 'Algorithm'))
    multiseq_dict['Singlet'] = multi_singlets
    multiseq_dict['Doublet'] = results['Doublet']
    multiseq_dict['Negative'] = results['Negative']
    multiseq_dict['Algorithm'] = results['Algorithm']
    return multiseq_dict

def htoDemul_visual(results):
    if results.get("Singlet") is None:   
        results['Singlet'] = 0
    return results

#This one gets a variable number of dicts from the algorithms used in the workflow
def create_df(*args):
    print(len(args))
    res_dict = defaultdict(list)
    for d in (args): # you can list as many input dicts as you want here
        for key, value in d.items():
            res_dict[key].append(value)
    return res_dict


def visualise(res_dict, path):
    df = pd.DataFrame(dict(res_dict))
    num_algorithms = len(res_dict)
    plot = df.plot.bar(x='Algorithm', stacked=True,color=["#DA4167","#E2B1B1","#258EA6"], title='Hashing Demultiplexing by algorithm',figsize=(15, 10))
    fig = plot.get_figure()
    graph = path+"Stacked_plot.png"
    fig.savefig(graph)



if __name__ == '__main__':
    hto = args.htoDemuxResultPath
    multi = args.multiSeqResultPath
    drops = args.hashedDropsResultPath

    hto_dict = pre_processing(hto,"htoDemux")
    hto_dict = htoDemul_visual(hto_dict)
    multi_pre = pre_processing(multi,"Multi-seq")
    multi_dict = multi_seq_visual(multi_pre)
    print(hto_dict)
    print(multi_dict)
    print("-----------")
    complete_dict = create_df(hto_dict,multi_dict)
    path = args.graphPath
    visualise(complete_dict,path)


