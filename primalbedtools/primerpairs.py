from primalbedtools.bedfiles import BedLine, StrandEnum


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

    def __init__(self, fbedlines: list[BedLine], rbedlines: list[BedLine], strict=True):
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

    def merge(self) -> tuple[BedLine, BedLine]:
        """
        Merges multiple forward and reverse primers into a single bedline
        """
        # Handle the fbedlines
        fbedline_start = min([bedline.start for bedline in self.fbedlines])
        fbedline_end = max([bedline.end for bedline in self.fbedlines])
        # Give it the longest seq to maintain the 7col format
        fbedline_sequence = max(
            [bedline.sequence for bedline in self.fbedlines], key=len
        )

        rbedline_start = min([bedline.start for bedline in self.rbedlines])
        rbedline_end = max([bedline.end for bedline in self.rbedlines])
        rbedline_sequence = max(
            [bedline.sequence for bedline in self.rbedlines], key=len
        )

        return BedLine(
            chrom=self.chrom,
            start=fbedline_start,
            end=fbedline_end,
            primername=f"{self.prefix}_{self.amplicon_number}_LEFT_1",
            pool=self.pool,
            strand=StrandEnum.FORWARD,
            sequence=fbedline_sequence,
        ), BedLine(
            chrom=self.chrom,
            start=rbedline_start,
            end=rbedline_end,
            primername=f"{self.prefix}_{self.amplicon_number}_RIGHT_1",
            pool=self.pool,
            strand=StrandEnum.REVERSE,
            sequence=rbedline_sequence,
        )
