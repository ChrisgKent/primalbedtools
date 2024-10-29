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
            print(
                f"All bedlines must have the same prefix ({','.join(prefixes)}). Using the first one."
            )
        self.prefix = sorted(prefixes)[0]

        # Check all chrom are the same
        chroms = set([bedline.chrom for bedline in all_lines])
        if len(chroms) != 1:
            raise ValueError(
                f"All bedlines must be on the same chromosome ({','.join(chroms)})"
            )
        self.chrom = chroms.pop()
        # Check all pools are the same
        pools = set([bedline.pool for bedline in all_lines])
        if len(pools) != 1:
            raise ValueError(
                f"All bedlines must be in the same pool ({','.join(map(str, pools))})"
            )
        self.pool = pools.pop()
        # Check all amplicon numbers are the same
        amplicon_numbers = set([bedline.amplicon_number for bedline in all_lines])
        if len(amplicon_numbers) != 1:
            raise ValueError(
                f"All bedlines must be the same amplicon ({','.join(map(str, amplicon_numbers))})"
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

    @property
    def amplicon_start(self) -> int:
        """Return the smallest start of the amplicon"""
        return min(self.fbedlines, key=lambda x: x.start).start

    @property
    def amplicon_end(self) -> int:
        """Return the largest end of the amplicon"""
        return max(self.rbedlines, key=lambda x: x.end).end

    @property
    def coverage_start(self) -> int:
        """Return the first base of coverage"""
        return max(self.fbedlines, key=lambda x: x.end).end

    @property
    def coverage_end(self) -> int:
        """Return the last base of coverage"""
        return min(self.rbedlines, key=lambda x: x.start).start

    @property
    def amplicon_name(self) -> str:
        """Return the name of the amplicon"""
        return f"{self.prefix}_{self.amplicon_number}"

    def to_amplicon_str(self) -> str:
        """Return the amplicon as a string in bed format"""
        return f"{self.chrom}\t{self.amplicon_start}\t{self.amplicon_end}\t{self.amplicon_name}\t{self.pool}"

    def to_primertrim_str(self) -> str:
        """Return the primertrimmed region as a string in bed format"""
        return f"{self.chrom}\t{self.coverage_start}\t{self.coverage_end}\t{self.amplicon_name}\t{self.pool}"


def create_primerpairs(bedlines: list[BedLine]) -> list[PrimerPair]:
    """
    Group bedlines into PrimerPair objects
    """
    grouped_bedlines = group_primer_pairs(bedlines)
    primer_pairs = []
    for fbedlines, rbedlines in grouped_bedlines:
        primer_pairs.append(PrimerPair(fbedlines, rbedlines))

    return primer_pairs
