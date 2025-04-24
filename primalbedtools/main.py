import argparse

from primalbedtools.bedfiles import (
    BedLineParser,
    merge_bedlines,
    sort_bedlines,
    update_primernames,
)
from primalbedtools.fasta import read_fasta
from primalbedtools.primerpairs import create_primerpairs
from primalbedtools.remap import remap
from primalbedtools.validate import validate, validate_primerbed


def main():
    parser = argparse.ArgumentParser(description="PrimalBedTools")
    subparsers = parser.add_subparsers(dest="subparser_name", required=True)

    # Remap subcommand
    remap_parser = subparsers.add_parser("remap", help="Remap BED file coordinates")
    remap_parser.add_argument("--bed", type=str, help="Input BED file", required=True)
    remap_parser.add_argument("--msa", type=str, help="Input MSA", required=True)
    remap_parser.add_argument(
        "--from_id", type=str, help="The ID to remap from", required=True
    )
    remap_parser.add_argument(
        "--to_id", type=str, help="The ID to remap to", required=True
    )

    # Sort subcommand
    sort_parser = subparsers.add_parser("sort", help="Sort BED file")
    sort_parser.add_argument("bed", type=str, help="Input BED file")

    # Update subcommand
    update_parser = subparsers.add_parser(
        "update", help="Update BED file with new information"
    )
    update_parser.add_argument("bed", type=str, help="Input BED file")

    # Amplicon subcommand
    amplicon_parser = subparsers.add_parser("amplicon", help="Create amplicon BED file")
    amplicon_parser.add_argument("bed", type=str, help="Input BED file")
    amplicon_parser.add_argument(
        "-t", "--primertrim", help="Primertrim the amplicons", action="store_true"
    )

    # merge subcommand
    merge_parser = subparsers.add_parser(
        "merge", help="Merge primer clouds into a single bedline"
    )
    merge_parser.add_argument("bed", type=str, help="Input BED file")

    # fasta subcommand
    fasta_parser = subparsers.add_parser("fasta", help="Convert .bed to .fasta")
    fasta_parser.add_argument("bed", type=str, help="Input BED file")

    # validate bedfile
    validate_bedfile_parser = subparsers.add_parser(
        "validate_bedfile", help="Validate a bedfile"
    )
    validate_bedfile_parser.add_argument("bed", type=str, help="Input BED file")

    # validate bedfile
    validate_parser = subparsers.add_parser(
        "validate", help="Validate a bedfile and reference"
    )
    validate_parser.add_argument("bed", type=str, help="Input BED file")
    validate_parser.add_argument("fasta", type=str, help="Input reference file")
    args = parser.parse_args()

    # Read in the bed file
    _headers, bedlines = BedLineParser.from_file(args.bed)

    if args.subparser_name == "remap":
        msa = read_fasta(args.msa)
        bedlines = remap(args.from_id, args.to_id, bedlines, msa)
    elif args.subparser_name == "sort":
        bedlines = sort_bedlines(bedlines)
    elif args.subparser_name == "update":
        bedlines = update_primernames(bedlines)
    elif args.subparser_name == "amplicon":
        primerpairs = create_primerpairs(bedlines)

        # Print the amplicons
        for primerpair in primerpairs:
            if args.primertrim:
                print(primerpair.to_primertrim_str())
            else:
                print(primerpair.to_amplicon_str())
        exit(0)  # Exit early
    elif args.subparser_name == "merge":
        bedlines = merge_bedlines(bedlines)
    elif args.subparser_name == "fasta":
        for line in bedlines:
            print(line.to_fasta(), end="")

        exit(0)  # Exit early
    elif args.subparser_name == "validate_bedfile":
        validate_primerbed(bedlines)
        exit(0)  # early exit

    elif args.subparser_name == "validate":
        validate(bedpath=args.bed, refpath=args.fasta)
        exit(0)  # early exit

    else:
        parser.print_help()

    bedfile_str = BedLineParser.to_str(_headers, bedlines)
    print(bedfile_str, end="")


if __name__ == "__main__":
    main()
