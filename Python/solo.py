import scvi
import scanpy as sc
import matplotlib.pyplot as plt
import anndata
import multiprocessing
from multiprocessing import freeze_support
import argparse

parser = argparse.ArgumentParser(description='Parser for SOLO - Doublet finding ')
#Input files
#H5 Available
parser.add_argument('--rna_data',  help='Input path to raw RNA expression matrix in 10x format.')
#Working good with 10x mtx
parser.add_argument('--soft',  help='Return probabilities instead of class label', default=False)
parser.add_argument('--max_epochs',  help='Max Epochs for training',type=int, default=100)
parser.add_argument('--lr',  help='Learning rate for training',type=float, default=0.001)

parser.add_argument('--output',  help='Output name',default="solo_prediction.csv")
args = parser.parse_args()

path_to_anndata = args.rna_data
adata = anndata.read_h5ad(path_to_anndata)
if __name__ == '__main__':
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train(max_epochs=args.max_epochs, lr=args.lr)
    prediction = solo.predict(soft=args.soft)
    prediction.to_csv(args.output)