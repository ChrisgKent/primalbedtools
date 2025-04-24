from primalbedtools.bedfiles import BedLine, BedLineParser
from primalbedtools.fasta import read_fasta
from primalbedtools.primerpairs import PrimerPair, create_primerpairs, do_pp_ol


def find_for_ol_in_pool(bedlines: list[BedLine]) -> set[tuple[PrimerPair, PrimerPair]]:
    """
    Sorts each bedline into its pool, and checks for overlap using the indexes.
    """
    # Create primerpairs
    primerpairs = create_primerpairs(bedlines)

    # keep track of each pp region
    regions = dict()
    # Get each
    for primerpair in primerpairs:
        # Add chrom
        if primerpair.chrom not in regions:
            regions[primerpair.chrom] = dict()

        # Add pool
        if primerpair.pool not in regions[primerpair.chrom]:
            regions[primerpair.chrom][primerpair.pool] = []

        # Add primer pair
        regions[primerpair.chrom][primerpair.pool].append(primerpair)

    ols = set()

    # For each check for ol
    for _chrom, pools in regions.items():
        for _pool, primerpairs in pools.items():
            for pp1 in primerpairs:
                for pp2 in primerpairs:
                    # Ignore self ol
                    if pp1 == pp2:
                        continue

                    # Add ordered ol to pool
                    if do_pp_ol(pp1, pp2):
                        d = tuple(sorted([pp1, pp2], key=lambda x: x.amplicon_start))
                        ols.add(d)  # type: ignore

    return ols


def validate_primerbed(bedlines: list[BedLine]):
    """
    This performs some simple QC on the structure and contents of a primer.bed file
    - Checks that each amplicon has at least one forward and reverse primer.
    - Checks for overlap between amplicons in the same pool.
    """
    # Check for overlaps in same pool
    # Also check for f and r primers
    overlaps = find_for_ol_in_pool(bedlines)
    if overlaps:
        ol_str = ", ".join(
            f"{pp1.amplicon_name}:{pp2.amplicon_name}" for pp1, pp2 in overlaps
        )
        raise ValueError(f"overlaps detected between: {ol_str}")


def validate_ref_and_bed(bedlines: list[BedLine], reference_path: str):
    fasta_ids = read_fasta(reference_path)

    # bedline chrom names
    bed_chrom_names = {bedline.chrom for bedline in bedlines}

    # Look for chroms in the bedfile that are not in the reference.fasta
    delta = bed_chrom_names - set(fasta_ids.keys())
    if delta:
        raise ValueError(f"chroms in primer.bed are not in reference.fasta: {delta}")
    # Look for chroms in the reference.fasta that are not in the bedfile
    delta = set(fasta_ids.keys()) - bed_chrom_names
    if delta:
        raise ValueError(f"chroms in reference.fasta are not in primer.bed: {delta}")


def validate(bedpath: str, refpath: str):
    # Read in bedlines
    ## Will validate bedline structure
    _header, bls = BedLineParser.from_file(bedpath)

    # Check ol, and left and right primer presence
    validate_primerbed(bls)

    # validate the bed and ref
    validate_ref_and_bed(bls, refpath)
