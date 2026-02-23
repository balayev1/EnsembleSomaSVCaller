#!/usr/bin/env python3

import argparse
import sys

def parse_args(args=None):
    Description = "Transform ascat purity/ploidy file to brass compatible input."
    Epilog = "Example usage: python transform_ascat_purityploidy.py <FILE_IN> <GENDER> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input ASCAT purityploidy.txt file")
    parser.add_argument("GENDER", help="Gender/Sex of the sample")
    parser.add_argument("FILE_OUT", help="Output formatted text file (male or female)")

    return parser.parse_args(args)

def make_outfile(file_in, gender, file_out):
    try:
        with open(file_in, 'r') as f:
            lines = [line.strip().split() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Input file {file_in} not found.", file=sys.stderr)
        sys.exit(1)

    if len(lines) < 2:
        print("Error: Invalid ASCAT purity file format.", file=sys.stderr)
        sys.exit(1)
    
    # Extract AberrantCellFraction (rho) and Ploidy
    rho = lines[1][0]
    ploidy = lines[1][1]

    if gender == 'female':
        gender_chr = "X"
        gender_found = "N"
    else:
        gender_chr = "Y"
        gender_found = "Y"

    # Write output
    with open(file_out, 'w') as out:
        out.write(f"rho {rho}\n")
        out.write(f"Ploidy {ploidy}\n")
        out.write(f"GenderChr {gender_chr}\n")
        out.write(f"GenderChrFound {gender_found}\n")

def main(args=None):
    args = parse_args(args)
    make_outfile(args.FILE_IN, args.GENDER, args.FILE_OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())