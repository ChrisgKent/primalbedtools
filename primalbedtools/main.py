import argparse

from primalbedtools.bedfiles import BedLineParser, sort_bedlines, update_primernames
from primalbedtools.fasta import read_fasta
from primalbedtools.remap import remap


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
    else:
        parser.print_help()

    bedfile_str = BedLineParser.to_str(_headers, bedlines)
    for line in bedfile_str.split("\n"):
        print(line)


if __name__ == "__main__":
    main()
