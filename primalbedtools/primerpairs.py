from primalbedtools.bedfiles import BedLine, group_primer_pairs


class PrimerPair:
    """
    A PrimerPair object represents an amplicon with forward and reverse primers.

    """

    fbedlines: list[BedLine]
    rbedlines: list[BedLine]

    chrom: str
    pool: int
    amplicon_number: int
    prefix: str

    def __init__(self, fbedlines: list[BedLine], rbedlines: list[BedLine]):
        self.fbedlines = fbedlines
        self.rbedlines = rbedlines

        all_lines = fbedlines + rbedlines

        # All prefixes must be the same
        prefixes = set([bedline.amplicon_prefix for bedline in all_lines])
        if len(prefixes) != 1:
            raise ValueError(
                f"All bedlines must have the same prefix, ({','.join(prefixes)})"
            )
        self.prefix = prefixes.pop()

        # Check all chrom are the same
        chroms = set([bedline.chrom for bedline in all_lines])
        if len(chroms) != 1:
            raise ValueError(
                f"All bedlines must be on the same chromosome, ({','.join(chroms)})"
            )
        self.chrom = chroms.pop()
        # Check all pools are the same
        pools = set([bedline.pool for bedline in all_lines])
        if len(pools) != 1:
            raise ValueError(
                f"All bedlines must be in the same pool, ({','.join(map(str, pools))})"
            )
        self.pool = pools.pop()
        # Check all amplicon numbers are the same
        amplicon_numbers = set([bedline.amplicon_number for bedline in all_lines])
        if len(amplicon_numbers) != 1:
            raise ValueError(
                f"All bedlines must be the same amplicon, ({','.join(map(str, amplicon_numbers))})"
            )
        self.amplicon_number = amplicon_numbers.pop()

        # Check both forward and reverse primers are present
        if not self.fbedlines:
            raise ValueError("No forward primers found")
        if not self.rbedlines:
            raise ValueError("No reverse primers found")

    @property
    def ipool(self) -> int:
        """Return the 0-based pool number"""
        return self.pool - 1

    @property
    def is_circular(self) -> bool:
        """Check if the amplicon is circular"""
        return self.fbedlines[0].end > self.fbedlines[0].start


def create_primerpairs(bedlines: list[BedLine]) -> list[PrimerPair]:
    """
    Group bedlines into PrimerPair objects
    """
    grouped_bedlines = group_primer_pairs(bedlines)
    primer_pairs = []
    for fbedlines, rbedlines in grouped_bedlines:
        primer_pairs.append(PrimerPair(fbedlines, rbedlines))

    return primer_pairs
