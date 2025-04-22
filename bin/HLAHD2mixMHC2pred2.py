#!/usr/bin/env python3

import sys
import re
import itertools
import argparse

def load_supported_alleles(alleles_file):
    """
    Loads supported HLA alleles from a file and creates a mapping between
    HLAHD format (e.g., HLA-DRB1*01:01) and MixMHC2pred format (e.g., DRB1_01_01).

    Args:
        alleles_file (str): Path to the file containing the list of supported alleles.
                            The file should have two tab-separated columns:
                            MixMHC2pred format and HLAHD format.

    Returns:
        dict: A dictionary mapping HLAHD format alleles to MixMHC2pred format alleles.
              Returns an empty dictionary if the file is empty or malformed.
    """
    supported_alleles = {}
    try:
        with open(alleles_file, "r", encoding="utf-8") as file:
            for line in file:
                parts = line.strip().split("\t")
                if len(parts) == 2 and "HLA-" in parts[1]:  # Ensure valid data
                    mixmhc2pred_format, official_name = parts
                    supported_alleles[official_name] = mixmhc2pred_format
    except FileNotFoundError:
        print(f"Error: Alleles file not found: {alleles_file}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"Error reading alleles file: {e}", file=sys.stderr)
        return {}
    return supported_alleles

def parse_hla_allele(allele):
    """
    Parses an HLA allele from HLAHD format (e.g., HLA-DRB1*01:01) and converts it
    to a simplified format (e.g., DRB1_01_01).

    Args:
        allele (str): The HLA allele in HLAHD format.

    Returns:
        str: The simplified allele format, or None if the allele is invalid or "Not typed".
    """
    if allele in ["Not typed", "-"]:
        return None
    match = re.match(r"HLA-([A-Z0-9]+)\*([\d:]+)", allele)
    if match:
        gene, allele_numbers = match.groups()
        allele_parts = allele_numbers.split(":")[:2]  # Keep only the first two fields
        return f"{gene}_{'_'.join(allele_parts)}"
    return None

def convert_hlahd_to_mixmhc2pred(input_file, output_file, alleles_file):
    """
    Converts HLAHD output to a format compatible with MixMHC2pred, using a list of
    supported alleles.

    Args:
        input_file (str): Path to the HLAHD output file.
        output_file (str): Path to the output file where the MixMHC2pred-compatible
                           alleles will be written.
        alleles_file (str): Path to the file containing the list of supported alleles.
    """
    # Load supported alleles
    supported_alleles = load_supported_alleles(alleles_file)
    if not supported_alleles:
        print("No supported alleles found. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Store parsed HLA class II alleles
    hla_class_ii = {"DRB1": [], "DRB3": [], "DRB4": [], "DRB5": [], "DQA1": [], "DQB1": [], "DPA1": [], "DPB1": []}

    try:
        # Read HLAHD output
        with open(input_file, 'r') as infile:
            for line in infile:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                gene, alleles = parts[0], parts[1:]

                if gene in hla_class_ii:
                    formatted_alleles = [parse_hla_allele(a) for a in alleles if parse_hla_allele(a)]
                    hla_class_ii[gene] = formatted_alleles
    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Generate DP and DQ combinations before filtering
    formatted_alleles = []

    # Add DRB1, DRB3, DRB4, DRB5 (listed individually)
    for drb_gene in ["DRB1", "DRB3", "DRB4", "DRB5"]:
        for allele in hla_class_ii[drb_gene]:
            if allele in supported_alleles.values():
                formatted_alleles.append(allele)

    # Generate all DQA1-DQB1 combinations
    dq_pairs = [
        f"{dqa}__{dqb}"
        for dqa, dqb in itertools.product(hla_class_ii["DQA1"], hla_class_ii["DQB1"])
    ]

    # Generate all DPA1-DPB1 combinations
    dp_pairs = [
        f"{dpa}__{dpb}"
        for dpa, dpb in itertools.product(hla_class_ii["DPA1"], hla_class_ii["DPB1"])
    ]

    # Filter DP and DQ combinations by supported alleles
    formatted_alleles.extend([pair for pair in dq_pairs if pair in supported_alleles.values()])
    formatted_alleles.extend([pair for pair in dp_pairs if pair in supported_alleles.values()])

    # Save output
    try:
        with open(output_file, 'w') as outfile:
            outfile.write("\n".join(formatted_alleles))
    except Exception as e:
        print(f"Error writing to output file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    Main function to parse command-line arguments and run the conversion process.
    """
    parser = argparse.ArgumentParser(
        description="Convert HLAHD output to MixMHC2pred-compatible format."
    )
    parser.add_argument("--input_file", help="Path to the HLAHD output file.")
    parser.add_argument("--output_file", help="Path to the output file.")
    parser.add_argument(
        "--alleles_file", help="Path to the Alleles_list_Human.txt file."
    )

    args = parser.parse_args()

    convert_hlahd_to_mixmhc2pred(
        args.input_file, args.output_file, args.alleles_file
    )
    print(f"Conversion complete. Output saved to {args.output_file}.")
    sys.exit(0)

if __name__ == "__main__":
    main()
