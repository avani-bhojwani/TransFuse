#!/usr/bin/env python

import pandas as pd
import argparse
import os
import glob

def merge_tpm(input_dir, output_file):
    merged_df = pd.DataFrame()
    quant_files = glob.glob(os.path.join(input_dir, '**/quant.sf'), recursive=True)

    for file_path in quant_files:
        # Extracting sample ID from the directory name
        sample_id = os.path.basename(os.path.dirname(file_path))

        # Read the quant.sf file
        df = pd.read_csv(file_path, sep='\t', usecols=['Name', 'TPM'])

        # Rename the 'TPM' column to the sample ID
        df.rename(columns={'TPM': sample_id}, inplace=True)

        # Merge with the main DataFrame
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='Name', how='outer')

    # Write the merged DataFrame to a file
    merged_df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Merge TPM values from Salmon quantification across samples")
    parser.add_argument("--input_dir", help="Directory containing sample folders with quant.sf files", required=True)
    parser.add_argument("-o", "--output_file", help="Output file name", required=True)
    
    args = parser.parse_args()
    merge_tpm(args.input_dir, args.output_file)

if __name__ == "__main__":
    main()