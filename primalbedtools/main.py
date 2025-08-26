import argparse
from importlib.metadata import version

from primalbedtools.amplicons import create_amplicons
from primalbedtools.bedfiles import (
    BedFileModifier,
)
from primalbedtools.fasta import read_fasta
from primalbedtools.remap import remap
from primalbedtools.scheme import Scheme
from primalbedtools.validate import validate, validate_primerbed


def main():
    parser = argparse.ArgumentParser(description="PrimalBedTools")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {version('primalbedtools')}"
    )

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

    # legacy parser
    downgrade_parser = subparsers.add_parser(
        "downgrade", help="Downgrade a bed file to an older version"
    )
    downgrade_parser.add_argument("bed", type=str, help="Input BED file")
    downgrade_parser.add_argument(
        "--merge-alts",
        help="Should alt primers be merged?",
        default=False,
        action="store_true",
    )
    # format
    format_parser = subparsers.add_parser("format", help="Format a bed file")
    format_parser.add_argument("bed", type=str, help="Input BED file")

    # format
    csv_parser = subparsers.add_parser("csv", help="Convert bed file to CSV")
    csv_parser.add_argument("bed", type=str, help="Input BED file")
    csv_parser.add_argument(
        "--no-headers", help="Remove the header row from the CSV", action="store_true"
    )
    csv_parser.add_argument(
        "--use-header-aliases",
        help="Should header aliases be used.",
        action="store_true",
    )

    args = parser.parse_args()

    # Read in the scheme
    scheme = Scheme.from_file(args.bed)

    if args.subparser_name == "remap":
        msa = read_fasta(args.msa)
        scheme.bedlines = remap(args.from_id, args.to_id, scheme.bedlines, msa)
        print(scheme.to_str(), end="")
        exit(0)
    elif args.subparser_name == "sort":
        scheme.bedlines = BedFileModifier.sort_bedlines(scheme.bedlines)
        print(scheme.to_str(), end="")
        exit(0)
    elif args.subparser_name == "update":
        scheme.bedlines = BedFileModifier.update_primernames(scheme.bedlines)
        print(scheme.to_str(), end="")
        exit(0)
    elif args.subparser_name == "amplicon":
        amplicons = create_amplicons(scheme.bedlines)

        # Print the amplicons
        for amplicon in amplicons:
            if args.primertrim:
                print(amplicon.to_primertrim_str())
            else:
                print(amplicon.to_amplicon_str())
        exit(0)  # Exit early
    elif args.subparser_name == "merge":
        scheme.bedlines = BedFileModifier.merge_primers(scheme.bedlines)
    elif args.subparser_name == "fasta":
        for line in scheme.bedlines:
            print(line.to_fasta(), end="")

        exit(0)  # Exit early
    elif args.subparser_name == "validate_bedfile":
        validate_primerbed(scheme.bedlines)
        exit(0)  # early exit

    elif args.subparser_name == "validate":
        validate(bedpath=args.bed, refpath=args.fasta)
        exit(0)  # early exit

    elif args.subparser_name == "downgrade":
        # merge primers if asked
        scheme.bedlines = BedFileModifier.downgrade_primernames(
            bedlines=scheme.bedlines, merge_alts=args.merge_alts
        )
        scheme.headers = []  # remove headers
        print(scheme.to_str(), end="")
        exit(0)

    elif args.subparser_name == "format":
        print(scheme.to_str(), end="")
        exit(0)
    elif args.subparser_name == "csv":
        print(
            scheme.to_delim_str(
                include_headers=not args.no_headers,
                use_header_aliases=args.use_header_aliases,
            )
        )
        exit(0)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
