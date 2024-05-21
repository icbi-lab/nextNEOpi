#!/usr/bin/env python

"""
Concatenate pvacseq outputfiles

Requirements:
    * Python >= 3.6.2
    * pandas >= 1.4

Copyright (c) 2024 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

import os
import argparse
import pandas as pd

def concat_files(directory, pattern):
    # List all files in the current directory that match the pattern
    files_to_concat = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(pattern)]
    print(files_to_concat)

    # Initialize an empty list to store DataFrames
    dfs = []

    # Iterate over the files and read them into DataFrames
    for file in files_to_concat:
        df = pd.read_csv(file, sep="\t")
        dfs.append(df)

    # Concatenate the DataFrames
    result = pd.concat(dfs, ignore_index=True)

    return result

if __name__ == "__main__":
    # Set up argparse to accept command-line arguments
    parser = argparse.ArgumentParser(description='Concatenate files in the current directory.')
    parser.add_argument('--dir', type=str, default='./', help='Directory to search for files (default: current directory)')
    parser.add_argument('--pattern', type=str, help='Pattern to match filenames')
    parser.add_argument('--output', type=str, default='out.tsv', help='Output filename')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to concatenate files with the specified pattern
    concatenated_df = concat_files(args.dir, args.pattern)

    # Save the concatenated DataFrame to a CSV file
    concatenated_df.to_csv(args.output, index=False, sep="\t", na_rep='NA')
