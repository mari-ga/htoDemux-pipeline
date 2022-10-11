from solo import hashsolo
import anndata
import matplotlib.pyplot as plt
import argparse
import pegasusio as io
import multiprocessing
from multiprocessing import freeze_support
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

parser = argparse.ArgumentParser(description='Parser for Hash Solo - Demultiplexing ')

parser.add_argument('--hto_data',  help='Input file filled only with hashing counts')
parser.add_argument('--priors_singlet',  help='Second element - prior for the singlet hypothesis', type=float,default=0.19)
parser.add_argument('--priors_doublet',  help='third element - prior for  the doublet hypothesis', type=float,default=0.8)
parser.add_argument('--priors_negative',  help='first element - prior for the negative hypothesis', type=float,default=0.01)

parser.add_argument('--number_of_noise_barcodes',  help='Number of threads to use.',type=int, default=None)
parser.add_argument('--output_file',  help='name for csv results',default="hash_solo.csv")
parser.add_argument('--output_plot',  help='name for the plot',default="hash_solo_plot.jpg")
args = parser.parse_args()

hto = args.hto_data
priors_list = list([args.priors_negative,args.priors_singlet,args.priors_doublet])

cell_hashing_data = io.read_input(hto)
if __name__ == '__main__':
    #multiprocessing.freeze_support()
    hashsolo.hashsolo(cell_hashing_data, priors=priors_list)
    cell_hashing_data.obs.to_csv(args.output_file)
    hashsolo.plot_qc_checks_cell_hashing(cell_hashing_data)
    plt.savefig(args.output_plot,dpi=400)






