from tokenize import String
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Parser for Summary Results')

parser.add_argument('--col_class_1',  help='name for the output file', default="--")
parser.add_argument('--col_2',  help='name for the output file', default="--")
parser.add_argument('--col_3',  help='name for the output file',default="--")
parser.add_argument('--col_4',  help='name for the output file',default="--")
parser.add_argument('--col_5',  help='name for the output file',default="--")
parser.add_argument('--col_6',  help='name for the output file',default="--")

parser.add_argument('--columns_assignment',  help='name for the output file',nargs='+', default=["--"])

parser.add_argument('--number_reports',  help='number of empty reports needed, always first classification name',default=1, type=int)

parser.add_argument('--output_report',  help='name for the output file',nargs='+', default=["report.csv"])

args = parser.parse_args()

columns = list()
columns.append(args.columns_classification)
columns.append(args.columns_assignment)


num_of_reports = args.number_reports
report_name = args.output_report

print(columns)
print(report_name)



for i in range(num_of_reports):
    temp_cols =  columns[i]
    df = pd.DataFrame(columns=temp_cols)
    temp_name = report_name[i]
    df.to_csv(temp_name)