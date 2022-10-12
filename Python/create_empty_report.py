from tokenize import String
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Parser for Summary Results')

parser.add_argument('--col_class_1',  help='name for the output file', default="--")
parser.add_argument('--col_2',  help='name for the output file', default="--")
parser.add_argument('--col_3',  help='name for the output file',default="--")



parser.add_argument('--output_report',  help='name for the output file',nargs='+', default=["empty.csv"])

args = parser.parse_args()

columns_class = list()
columns_assign = list()
columns_class.append(args.columns_classification)
columns_assign.append(args.columns_assignment)


num_of_reports = args.number_reports
report_name = args.output_report

if num_of_reports==1:
    
    print("ok")

else:
    print ("HTODemux")
    columns_assign = ["Barcode_htoDemux","Assignment-htoDemux"]
    columns_class = ["Barcode_htoDemux","Classification-HTODemux"]
    report_class = pd.DataFrame(columns=columns_class)
    report_assign = pd.DataFrame(columns=columns_assign)
    report_class.to_csv('')
    report_assign.to_csv('')

for i in range(num_of_reports):
    temp_cols =  columns[i]
    df = pd.DataFrame(columns=temp_cols)
    temp_name = report_name[i]
    df.to_csv(temp_name)