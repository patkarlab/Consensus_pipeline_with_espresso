#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Take merged output (varscan, flt3r, coverage) and get_itd output, find common variants between flt3r and get_itd, add to combined Excel file")

parser.add_argument("--getitd_out", required=True, help="Path to get_itd output")
parser.add_argument("--merged_out", required=True, help="Path to merged Excel file (contains varscan, filt3r, coverage sheets)")
parser.add_argument("--output", required=True, help="Path to output Excel file")

args = parser.parse_args()

# read in input files from merged Excel
merged_sheet1 = pd.read_excel(args.merged_out, sheet_name="varscan")
df1 = pd.read_excel(args.merged_out, sheet_name="filt3r")
cov_bed = pd.read_excel(args.merged_out, sheet_name="coverage")
df2 = pd.read_csv(args.getitd_out, sep="\t")

# add tool name suffix to each column
df1 = df1.add_suffix('_flt3r')
df2 = df2.add_suffix('_getitd')

# join the 2 dataframes on the basis of key, drop the key column
df3 = df1.assign(key=0).merge(df2.assign(key=0), how='left', on='key').drop(columns=['key'])

# function to find the longest common substring(min length=5) between 2 sequences
def common_substring(seq1, seq2):
    min_len = 5
    max_len = 0
    match_seq = ""

    if isinstance(seq1, str) and isinstance(seq2, str) and len(seq1) == len(seq2):
        for i in range(len(seq1)):
            for j in range(1 + min_len, len(seq1) + 1):
                substr = seq1[i:j]
                if substr in seq2 and len(substr) >= min_len:
                    if len(substr) > max_len:
                        max_len = len(substr)
                        match_seq = substr

    return match_seq if len(match_seq) > min_len else False

try:
    required_cols = ['sequence_flt3r', 'seq_getitd']

    if not df3.empty:
        if all(col in df3.columns for col in required_cols) and \
           not df3[required_cols].isna().all().any():

            df3['matching_itd'] = df3.apply(
                lambda row: common_substring(row["sequence_flt3r"], row["seq_getitd"]), axis=1
            )
        else:
            print("Skipping processing: Required columns exist but contain only NaNs.")
            df3['matching_itd'] = False
    else:
        print("Skipping processing: Merged DataFrame is completely empty.")
        df3['matching_itd'] = False

except Exception as e:
    print(f"Error during processing: {e}")
    if 'matching_itd' not in df3.columns:
        df3['matching_itd'] = False

# move the common sequence column to index 0
common_itd = df3.pop('matching_itd')
df3.insert(0, 'matching_itd', common_itd)

# extract rows where common variant identified
match_df = df3[df3['matching_itd'] != False]

# write output
with pd.ExcelWriter(args.output) as writer:
    merged_sheet1.to_excel(writer, sheet_name="varscan", index=False)
    df1.to_excel(writer, sheet_name="flt3r", index=False)
    df2.to_excel(writer, sheet_name="get_itd", index=False)
    match_df.to_excel(writer, sheet_name="common_itds", index=False)
    cov_bed.to_excel(writer, sheet_name="coverage", index=False)
