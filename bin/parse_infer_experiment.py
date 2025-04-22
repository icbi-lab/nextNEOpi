#!/usr/bin/env python3
import sys
import re
import argparse

def parse_infer_experiment(file_path):
    """
    Parses infer_experiment.py output to determine the strandedness.

    Args:
        file_path (str): Path to the infer_experiment.py output file.

    Returns:
        str: "RF" (Reverse stranded), "FR" (Forward stranded), or "XS" (unstranded).
             Returns None if the file is not found or if parsing fails.
    """
    try:
        fr_fraction, rf_fraction = None, None

        with open(file_path, "r") as file:
            for line in file:
                # Use more robust regex patterns to handle potential variations in output
                fr_match = re.search(r"Fraction of reads explained by \"1\+\+,1--\,2\+\-\,2\-\+\"\:\s+([\d\.]+)", line)
                rf_match = re.search(r"Fraction of reads explained by \"1\+\-\,1\-\+\,2\+\+,2--\"\:\s+([\d\.]+)", line)

                if fr_match:
                    fr_fraction = float(fr_match.group(1))
                elif rf_match:
                    rf_fraction = float(rf_match.group(1))

        if fr_fraction is None or rf_fraction is None:
            raise ValueError("Error parsing file: Missing required fractions.")

        # Determine strandedness based on thresholds
        if fr_fraction >= 0.8:
            return "FR"  # Forward stranded
        elif rf_fraction >= 0.8:
            return "RF"  # Reverse stranded
        else:
            return "XS"  # Unstranded

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}", file=sys.stderr)
        return None
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        return None


def main():
    """
    Main function to handle command-line arguments and execute the parsing.
    """
    parser = argparse.ArgumentParser(
        description="Parse infer_experiment.py output to determine strandedness."
    )
    parser.add_argument(
        "--infer_experiment_output",
        dest="file_path",
        required=True,
        metavar="FILE_PATH",
        help="Path to the infer_experiment.py output file.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output (print fractions).",
    )

    args = parser.parse_args()

    result = parse_infer_experiment(args.file_path)

    if result:
        print(result)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
