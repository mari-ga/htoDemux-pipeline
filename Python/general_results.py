import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from matplotlib import cm
from upsetplot import UpSet
import argparse

parser = argparse.ArgumentParser(description='Parser for General results plotting')

parser.add_argument('--htodemul_results',  help='CSV file with classification results from HTOdemul')
parser.add_argument('--htodemul_assignment',  help='CSV file with Hashtag assignment results from HTOdemul')
parser.add_argument('--multiseq_results',  help='CSV file with classification results from MULTI-seq')
parser.add_argument('--hashed_drop_results',  help='CSV file with classification results from Hashed Drop')
parser.add_argument('--demuxem_results',  help='CSV file with classification results from DemuxEM')
parser.add_argument('--hashsolo_results',  help='CSV file with classification results from MULTI-seq')

args = parser.parse_args()


