import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Parser for Summary Results')

parser.add_argument('--report_name',  help='Report name', default="report.csv")

args = parser.parse_args()


total_classification = pd.DataFrame(list())
total_classification.to_csv(args.report_name)