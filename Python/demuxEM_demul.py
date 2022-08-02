from ast import arg
import string
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pegasusio as io
from demuxEM.tools import *
import pegasus
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser(description='Parser for DemuxEM - Demultiplexing ')
#Input files
parser.add_argument('--rna_data',  help='Input raw RNA expression matrix in 10x hdf5 format.')
parser.add_argument('--hto_matrix',  help='HTO (antibody tag) count matrix in CSV format.')
#output name
#parser.add_argument('--output',  help='Output name',default="demuxEm")
#Parameters
parser.add_argument('--threads',  help='Number of threads',default=1)
parser.add_argument('--alpha',  help='The Dirichlet prior concentration parameter (alpha) on samples.', default=0.0)
parser.add_argument('--alpha_noise',  help='The Dirichlet prior concentration parameter (alpha) on samples.', default=1.0)
parser.add_argument('--min_signal',  help='Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown.', default=10.0)
parser.add_argument('--tol',  help='Threshold used for the EM convergence..', default=1e-6)
#Files - plots

args = parser.parse_args()


rna_matrix = args.rna_data
hto_matrix = args.hto_matrix

pegasus.demultiplex(rna_matrix,hto_matrix,min_signal=args.min_signal, alpha=args.alpha,alpha_noise=args.alpha_noise,tol=args.tol,n_threads= args.threads )

