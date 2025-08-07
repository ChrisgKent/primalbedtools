from typing import Optional

from primalbedtools.bedfiles import BedLine, PrimerClass, group_amplicons


class Amplicon:
    """
    A Amplicon object represents an PCR amplicon with forward and reverse primers, and optional probes

    """

    left: list[BedLine]
    right: list[BedLine]
    probes: list[BedLine]

    chrom: str
    pool: int
    amplicon_number: int
    prefix: str

    def __init__(
        self,
        left: list[BedLine],
        right: list[BedLine],
        probes: Optional[list[BedLine]] = None,
    ):
        self.left = left
        self.right = right

        if probes is None:
            probes = []
        self.probes = probes

        all_lines = left + right + probes

        # All prefixes must be the same
        prefixes = set([bedline.amplicon_prefix for bedline in all_lines])
        prefixes = sorted(prefixes)

        if len(prefixes) != 1:
            print(
                f"All bedlines must have the same prefix ({','.join(prefixes)}). Using the alphanumerically first one ({prefixes[0]})."
            )
        self.prefix = prefixes[0]

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
        if not self.left:
            raise ValueError(
                f"No forward primers found for {self.prefix}_{self.amplicon_number}"
            )
        if not self.right:
            raise ValueError(
                f"No reverse primers found for {self.prefix}_{self.amplicon_number}"
            )

    @property
    def ipool(self) -> int:
        """Return the 0-based pool number"""
        return self.pool - 1

    @property
    def is_circular(self) -> bool:
        """Check if the amplicon is circular"""
        return self.left[0].end > self.right[0].start

    @property
    def amplicon_start(self) -> int:
        """Return the smallest start of the amplicon"""
        return min(self.left, key=lambda x: x.start).start

    @property
    def amplicon_end(self) -> int:
        """Return the largest end of the amplicon"""
        return max(self.right, key=lambda x: x.end).end

    @property
    def coverage_start(self) -> int:
        """Return the first base of coverage"""
        return max(self.left, key=lambda x: x.end).end

    @property
    def coverage_end(self) -> int:
        """Return the last base of coverage"""
        return min(self.right, key=lambda x: x.start).start

    @property
    def amplicon_name(self) -> str:
        """Return the name of the amplicon"""
        return f"{self.prefix}_{self.amplicon_number}"

    @property
    def probe_region(self) -> Optional[tuple[int, int]]:
        """
        Returns the half open position of the PROBES (if present).
        """
        if not self.probes:
            return None
        return (min(p.start for p in self.probes), max(p.end for p in self.probes))

    @property
    def left_region(self) -> tuple[int, int]:
        """
        Returns the half open position of the LEFT primers
        """
        return (min(lp.start for lp in self.left), max(lp.end for lp in self.left))

    @property
    def right_region(self) -> tuple[int, int]:
        """
        Returns the half open position of the RIGHT primers
        """
        return (min(rp.start for rp in self.right), max(rp.end for rp in self.right))

    def to_amplicon_str(self) -> str:
        """Return the amplicon as a string in bed format"""
        return f"{self.chrom}\t{self.amplicon_start}\t{self.amplicon_end}\t{self.amplicon_name}\t{self.pool}"

    def to_primertrim_str(self) -> str:
        """Return the primertrimmed region as a string in bed format"""
        return f"{self.chrom}\t{self.coverage_start}\t{self.coverage_end}\t{self.amplicon_name}\t{self.pool}"


def create_amplicons(bedlines: list[BedLine]) -> list[Amplicon]:
    """
    Group bedlines into Amplicon objects
    """
    grouped_bedlines = group_amplicons(bedlines)
    primer_pairs = []
    for pdict in grouped_bedlines:
        primer_pairs.append(
            Amplicon(
                left=pdict.get(PrimerClass.LEFT.value, []),
                right=pdict.get(PrimerClass.RIGHT.value, []),
                probes=pdict.get(PrimerClass.PROBE.value, []),
            )
        )

    return primer_pairs


def do_pp_ol(pp1: Amplicon, pp2: Amplicon) -> bool:
    if range(
        max(pp1.amplicon_start, pp2.amplicon_start),
        min(pp1.amplicon_end, pp2.amplicon_end) + 1,
    ):
        return True
    else:
        return False
