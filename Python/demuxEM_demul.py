import pegasusio as io
from demuxEM.tools import *
import pegasus as pg
import argparse
import multiprocessing
from multiprocessing import freeze_support

parser = argparse.ArgumentParser(description='Parser for DemuxEM - Demultiplexing ')
#Input files
#H5 Available
parser.add_argument('--rna_data',  help='Input raw RNA expression matrix in 10x hdf5 format.')
#Working good with 10x mtx
parser.add_argument('--hto_matrix',  help='HTO (antibody tag) count matrix in mtx format.')
#output name
parser.add_argument('--alpha',  help='The Dirichlet prior concentration parameter (alpha) on samples.', type=float, default=0.0)
parser.add_argument('--alpha_noise',  help='The Dirichlet prior concenration parameter on the background noise',type=float, default=1.0)
parser.add_argument('--tol',  help='Threshold used for the EM convergence',type=float, default=1e-6)
parser.add_argument('--n_threads',  help='Number of threads to use.',type=int, default=1)
parser.add_argument('--min_signal',  help='Any cell/nucleus with less than min_signal hashtags from the signal will be marked as Negative.',type=float, default=10.0)

parser.add_argument('--output',  help='Output name',default="demuxEm.csv")

#Files - plots

args = parser.parse_args()

umi_mat = args.rna_data
hto_mat = args.hto_matrix
alpha_val = args.alpha
alpha_noise_val = args.alpha_noise
min_signal_val = args.min_signal
tol_val = args.tol
n_threads_val = args.n_threads


#The input is read with 2 libraries, only pegasus and pegasus io for the mtx
#Anyways, both matrices are transformed to multimodal
hto_data = io.read_input(hto_mat)
umi_data = pg.read_input(input_file=umi_mat)
if __name__ == '__main__':
    multiprocessing.freeze_support()
    pg.qc_metrics(umi_data)
    pg.identify_robust_genes(umi_data)
    pg.estimate_background_probs(hto_data)
    pg.demultiplex(umi_data, hto_data,alpha=alpha_val, alpha_noise=alpha_noise_val,min_signal=min_signal_val,tol=tol_val,n_threads=n_threads_val)
    umi_data.obs.to_csv(args.output)
    

