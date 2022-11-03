import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from matplotlib import cm
import argparse



parser = argparse.ArgumentParser(description='Parser for Final Report')

parser.add_argument('--general_classification',  help='CSV file with General Classification results', default="-")
parser.add_argument('--solo_results',  help='CSV file with prediction results from Solo',default="-")
parser.add_argument('--output_final',  help='name for the output file',default="final_report.csv")

args = parser.parse_args()



total_classification = pd.DataFrame()

if(args.general_classification != "-" ):
    general = pd.read_csv(args.general_classification)
    total_classification = pd.concat([total_classification,general],axis=1)



if(args.solo_results != "-" ):
    solo_results = pd.read_csv(args.solo_results)
    if solo_results.shape[1] == 2:
            solo_results.columns = ["Barcode_Solo","Prediction-Solo"]
            total_classification = pd.concat([total_classification,solo_results],axis=1)
    else:
            solo_results.columns = ["Barcode_Solo","Doublet", "Singlet"]
            total_classification = pd.concat([total_classification,solo_results],axis=1)


total_classification.to_csv(args.output_final)   