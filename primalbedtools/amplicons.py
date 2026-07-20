import logging
from typing import Optional

from primalbedtools.bedfiles import (
    BedLine,
    BedLineParser,
    PrimerClass,
    group_amplicons,
)
from primalbedtools.logs import get_logger

logger = get_logger(__name__)

# Cap on how many amplicon names a divergence summary spells out per prefix set
MAX_REPORTED_AMPLICONS = 10


def _group_by_value(bedlines: list[BedLine], attr: str, label: str) -> str:
    """Render bedline primernames grouped by the value of ``attr``, one line per value.

    Produces a human-readable summary such as::

        pool 1: scheme_1_LEFT_1, scheme_1_RIGHT_1
        pool 2: scheme_1_LEFT_2

    so error messages can name exactly which primers diverge.
    """
    groups: dict = {}
    for bl in bedlines:
        groups.setdefault(getattr(bl, attr), []).append(bl.primername)
    return "\n".join(
        f"  {label} {value}: {', '.join(sorted(names))}"
        for value, names in sorted(groups.items(), key=lambda kv: str(kv[0]))
    )


def _format_prefix_divergence(
    amplicons: list["Amplicon"],
    max_reported: int = MAX_REPORTED_AMPLICONS,
) -> str:
    """Summarise amplicons whose primers carry more than one prefix.

    Amplicons are grouped by their exact set of prefixes, so a scheme-wide split
    (every amplicon mixing "<scheme>-F" with "<scheme>-R", say) collapses to a
    single line rather than one line per amplicon.
    """
    groups: dict = {}
    for amplicon in amplicons:
        groups.setdefault(tuple(amplicon.prefixes), []).append(amplicon)

    lines = [
        f"{len(amplicons)} amplicon(s) have primers with more than one prefix. "
        "The joined prefix is used in the amplicon name."
    ]
    for prefixes, group in sorted(groups.items()):
        group = sorted(group, key=lambda a: (a.chrom, a.amplicon_number))
        shown = ", ".join(a.amplicon_name for a in group[:max_reported])
        remaining = len(group) - max_reported
        if remaining > 0:
            shown += f", ... (+{remaining} more)"
        lines.append(f"  prefixes {', '.join(prefixes)}: {shown}")
    return "\n".join(lines)


class Amplicon:
    """A class representing a PCR amplicon with forward and reverse primers, and optional probes.

    An Amplicon object encapsulates all primers (LEFT, RIGHT, and optional PROBE) that belong
    to the same amplicon number and pool. It provides methods to calculate amplicon boundaries,
    coverage regions, and validate that all primers are consistent.

    Attributes:
        left (list[BedLine]): List of LEFT primer BedLine objects
        right (list[BedLine]): List of RIGHT primer BedLine objects
        probes (list[BedLine]): List of PROBE primer BedLine objects
        chrom (str): Chromosome name where the amplicon is located
        pool (int): 1-based pool number
        amplicon_number (int): Amplicon number from primer names
        prefix (str): Amplicon prefix from primer names. If the primers carry
            more than one prefix, all of them joined with "-"
        prefixes (list[str]): Sorted, de-duplicated prefixes found on the
            primers. More than one entry means the primers disagree

    Raises:
        ValueError: If primers have inconsistent chromosome, pool, or amplicon numbers
        ValueError: If LEFT or RIGHT primers are missing

    Examples:
        >>> from primalbedtools.bedfiles import BedLine
        >>> left_primer = BedLine(chrom="chr1", start=100, end=120,
        ...                       primername="scheme_1_LEFT_1", pool=1,
        ...                       strand="+", sequence="ACGT")
        >>> right_primer = BedLine(chrom="chr1", start=200, end=220,
        ...                        primername="scheme_1_RIGHT_1", pool=1,
        ...                        strand="-", sequence="ACGT")
        >>> amplicon = Amplicon(left=[left_primer], right=[right_primer])
        >>> print(amplicon.amplicon_name)
        scheme_1
        >>> print(amplicon.amplicon_start, amplicon.amplicon_end)
        100 220
    """

    left: list[BedLine]
    right: list[BedLine]
    probes: list[BedLine]

    chrom: str
    pool: int
    amplicon_number: int
    prefix: str
    prefixes: list[str]

    def __init__(
        self,
        left: list[BedLine],
        right: list[BedLine],
        probes: Optional[list[BedLine]] = None,
    ):
        """Initialize an Amplicon with LEFT and RIGHT primers, and optional PROBE primers.

        Args:
            left: List of BedLine objects representing LEFT primers
            right: List of BedLine objects representing RIGHT primers
            probes: Optional list of BedLine objects representing PROBE primers

        Raises:
            ValueError: If primers have inconsistent chromosome, pool, or amplicon numbers
            ValueError: If LEFT or RIGHT primers are missing
        """
        self.left = left
        self.right = right

        if probes is None:
            probes = []
        self.probes = probes

        all_lines = left + right + probes

        # Prefixes should agree. When they don't, every prefix is kept in the
        # name so the divergence stays visible rather than being silently
        # resolved to one of them.
        prefixes = sorted({bedline.amplicon_prefix for bedline in all_lines})

        if len(prefixes) > 1 and logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                f"Amplicon primers have differing prefixes ({', '.join(prefixes)}); "
                f"using the joined prefix ({'-'.join(prefixes)}).\n"
                + _group_by_value(all_lines, "amplicon_prefix", "prefix")
            )

        self.prefixes = prefixes
        self.prefix = "-".join(prefixes)

        # Check all chrom are the same
        chroms = set([bedline.chrom for bedline in all_lines])
        if len(chroms) != 1:
            raise ValueError(
                "Failed to create amplicon as provided bedlines are on different chromosomes:\n"
                + _group_by_value(all_lines, "chrom", "chrom")
                + "\n\n"
                + BedLineParser.to_str(bedlines=all_lines, headers=None)
            )
        self.chrom = chroms.pop()

        # Check all amplicon numbers are the same
        amplicon_numbers = set([bedline.amplicon_number for bedline in all_lines])
        if len(amplicon_numbers) != 1:
            raise ValueError(
                "Failed to create amplicon as provided bedlines have different amplicon numbers:\n"
                + _group_by_value(all_lines, "amplicon_number", "amplicon")
                + "\n\n"
                + BedLineParser.to_str(bedlines=all_lines, headers=None)
            )
        self.amplicon_number = amplicon_numbers.pop()

        # Check all pools are the same
        pools = set([bedline.pool for bedline in all_lines])
        if len(pools) != 1:
            raise ValueError(
                f"Failed to create amplicon ({self.prefix}_{self.amplicon_number}) as provided bedlines have different pools:\n"
                + _group_by_value(all_lines, "pool", "pool")
                + "\n\n"
                + BedLineParser.to_str(bedlines=all_lines, headers=None)
            )
        self.pool = pools.pop()

        # Check both forward and reverse primers are present
        if not self.left:
            present = ", ".join(sorted(bl.primername for bl in self.right)) or "none"
            raise ValueError(
                f"Failed to create amplicon ({self.prefix}_{self.amplicon_number}) as no forward primers found. "
                f"Reverse primers present: {present}\n\n"
                + BedLineParser.to_str(bedlines=all_lines, headers=None)
            )
        if not self.right:
            present = ", ".join(sorted(bl.primername for bl in self.left)) or "none"
            raise ValueError(
                f"Failed to create amplicon ({self.prefix}_{self.amplicon_number}) as no reverse primers found. "
                f"Forward primers present: {present}\n\n"
                + BedLineParser.to_str(bedlines=all_lines, headers=None)
            )

    def __lt__(self, other):
        if not isinstance(other, Amplicon):
            return NotImplemented
        return (self.chrom, self.amplicon_number, self.pool) < (
            other.chrom,
            other.amplicon_number,
            other.pool,
        )

    def __le__(self, other):
        if not isinstance(other, Amplicon):
            return NotImplemented
        return (self.chrom, self.amplicon_number, self.pool) <= (
            other.chrom,
            other.amplicon_number,
            other.pool,
        )

    def __gt__(self, other):
        if not isinstance(other, Amplicon):
            return NotImplemented
        return (self.chrom, self.amplicon_number, self.pool) > (
            other.chrom,
            other.amplicon_number,
            other.pool,
        )

    def __ge__(self, other):
        if not isinstance(other, Amplicon):
            return NotImplemented
        return (self.chrom, self.amplicon_number, self.pool) >= (
            other.chrom,
            other.amplicon_number,
            other.pool,
        )

    def __eq__(self, other):
        if not isinstance(other, Amplicon):
            return NotImplemented
        return (self.chrom, self.amplicon_number, self.pool) == (
            other.chrom,
            other.amplicon_number,
            other.pool,
        )

    def __ne__(self, other):
        if not isinstance(other, Amplicon):
            return NotImplemented
        return not self.__eq__(other)

    def __hash__(self) -> int:
        """Hash is based off the self.to_amplicon_str()"""
        return hash(self.to_amplicon_str())

    @property
    def ipool(self) -> int:
        """Get the 0-based pool number.

        Returns:
            int: Pool number converted from 1-based to 0-based indexing
        """
        return self.pool - 1

    @property
    def is_circular(self) -> bool:
        """Check if the amplicon appears to be circular (LEFT primer end > RIGHT primer start).

        This can indicate a circular genome where the amplicon spans the origin.

        Returns:
            bool: True if any LEFT primer end is greater than any RIGHT primer start
        """
        return self.left[0].end > self.right[0].start

    @property
    def amplicon_start(self) -> int:
        """Get the start position of the amplicon (earliest LEFT primer start).

        Returns:
            int: The smallest start position among all LEFT primers
        """
        return min(self.left, key=lambda x: x.start).start

    @property
    def amplicon_end(self) -> int:
        """Get the end position of the amplicon (latest RIGHT primer end).

        Returns:
            int: The largest end position among all RIGHT primers
        """
        return max(self.right, key=lambda x: x.end).end

    @property
    def coverage_start(self) -> int:
        """Get the start of the coverage region (latest LEFT primer end).

        This represents the first base that would be covered after primer trimming.

        Returns:
            int: The largest end position among all LEFT primers
        """
        return max(self.left, key=lambda x: x.end).end

    @property
    def coverage_end(self) -> int:
        """Get the end of the coverage region (earliest RIGHT primer start).

        This represents the last base that would be covered after primer trimming.

        Returns:
            int: The smallest start position among all RIGHT primers
        """
        return min(self.right, key=lambda x: x.start).start

    @property
    def amplicon_name(self) -> str:
        """Get the name of the amplicon.

        Returns:
            str: Amplicon name in format "prefix_number"
        """
        return f"{self.prefix}_{self.amplicon_number}"

    @property
    def probe_region(self) -> Optional[tuple[int, int]]:
        """Get the genomic region covered by PROBE primers.

        Returns:
            Optional[tuple[int, int]]: Half-open interval (start, end) of PROBE region,
                                     or None if no probes are present
        """
        if not self.probes:
            return None
        return (
            min(p.start for p in self.probes),
            max(p.end for p in self.probes),
        )

    @property
    def left_region(self) -> tuple[int, int]:
        """Get the genomic region covered by LEFT primers.

        Returns:
            tuple[int, int]: Half-open interval (start, end) of LEFT primer region
        """
        return (
            min(lp.start for lp in self.left),
            max(lp.end for lp in self.left),
        )

    @property
    def right_region(self) -> tuple[int, int]:
        """Get the genomic region covered by RIGHT primers.

        Returns:
            tuple[int, int]: Half-open interval (start, end) of RIGHT primer region
        """
        return (
            min(rp.start for rp in self.right),
            max(rp.end for rp in self.right),
        )

    def to_amplicon_str(self) -> str:
        """Convert the amplicon to a BED format string representing the full amplicon.

        Returns:
            str: Tab-delimited string with chrom, amplicon_start, amplicon_end,
                 amplicon_name, and pool
        """
        return f"{self.chrom}\t{self.amplicon_start}\t{self.amplicon_end}\t{self.amplicon_name}\t{self.pool}"

    def to_primertrim_str(self) -> str:
        """Convert the amplicon to a BED format string representing the coverage region.

        This represents the region that would remain after primer trimming.

        Returns:
            str: Tab-delimited string with chrom, coverage_start, coverage_end,
                 amplicon_name, and pool
        """
        return f"{self.chrom}\t{self.coverage_start}\t{self.coverage_end}\t{self.amplicon_name}\t{self.pool}"


def create_amplicons(bedlines: list[BedLine]) -> list[Amplicon]:
    """Group BedLine objects into Amplicon objects by chromosome, amplicon number, and pool.

    Args:
        bedlines: List of BedLine objects to group into amplicons

    Returns:
        list[Amplicon]: List of Amplicon objects created from the input bedlines

    Raises:
        ValueError: If any amplicon is missing LEFT or RIGHT primers
        ValueError: If primers within an amplicon have inconsistent attributes

    Note:
        Primers that disagree on their prefix are not an error. A single warning
        summarising every affected amplicon is logged, and each such amplicon is
        named using all of its prefixes joined with "-".

    Examples:
        >>> from primalbedtools.bedfiles import BedLineParser
        >>> headers, bedlines = BedLineParser.from_file("primers.bed")
        >>> amplicons = create_amplicons(bedlines)
        >>> print(f"Created {len(amplicons)} amplicons")
    """
    grouped_bedlines = group_amplicons(bedlines)
    primer_pairs = []
    divergent = []
    for pdict in grouped_bedlines:
        amplicon = Amplicon(
            left=pdict.get(PrimerClass.LEFT.value, []),
            right=pdict.get(PrimerClass.RIGHT.value, []),
            probes=pdict.get(PrimerClass.PROBE.value, []),
        )
        primer_pairs.append(amplicon)
        if len(amplicon.prefixes) > 1:
            divergent.append(amplicon)

    # Reported once for the whole scheme rather than once per amplicon
    if divergent:
        logger.warning(_format_prefix_divergence(divergent))

    return primer_pairs


def do_pp_ol(pp1: Amplicon, pp2: Amplicon) -> bool:
    """Check if two amplicons have overlapping genomic regions.

    Args:
        pp1: First amplicon to compare
        pp2: Second amplicon to compare

    Returns:
        bool: True if the amplicons have any overlapping genomic coordinates

    Examples:
        >>> amplicon1 = Amplicon(left=[...], right=[...])  # Covers 100-200
        >>> amplicon2 = Amplicon(left=[...], right=[...])  # Covers 150-250
        >>> do_pp_ol(amplicon1, amplicon2)
        True
    """
    return max(pp1.amplicon_start, pp2.amplicon_start) < min(
        pp1.amplicon_end, pp2.amplicon_end
    )
