from difflib import ndiff, unified_diff
from typing import Optional

from primalbedtools.bedfiles import BedLine, BedLineParser

PLACEHOLDER_STR = "<placeholder>"


def create_normalised_bedfile_str(
    bedlines1: list[BedLine],
    bedlines2: list[BedLine],
    header1: Optional[list[str]] = None,
    header2: Optional[list[str]] = None,
    ignore_order: bool = True,
    ignore_attr: bool = False,
    ignore_header: bool = False,
    ignore_primer_prefix: bool = False,
) -> tuple[str, str]:
    """Creates normalised bedfile strings formatted for diff comparison.

    Generates string representations of two sets of bedlines, optionally normalizing
    order, attributes, and headers to facilitate meaningful comparisons.

    Args:
        bedlines1 (list[BedLine]): The first list of BedLine objects.
        bedlines2 (list[BedLine]): The second list of BedLine objects.
        header1 (Optional[list[str]] , optional): Headers for the first bedfile. Defaults to None.
        header2 (Optional[list[str]] , optional): Headers for the second bedfile. Defaults to None.
        ignore_order (bool, optional): If True, sorts bedlines before string generation. Defaults to True.
        ignore_attr (bool, optional): If True, excludes attributes from the string representation. Defaults to False.
        ignore_header (bool, optional): If True, excludes headers from the string representation. Defaults to False.
        ignore_primer_prefix (bool, optional): If True, excludes the primername prefix from the string representation. Defaults to False.


    Returns:
        tuple[str, str]: A tuple containing the two normalised bedfile strings.
    """
    if ignore_header:
        header1 = None
        header2 = None

    if ignore_order:
        bed_str1 = BedLineParser.to_str(header1, sorted(bedlines1), ignore_attr)
        bed_str2 = BedLineParser.to_str(header2, sorted(bedlines2), ignore_attr)
    else:
        bed_str1 = BedLineParser.to_str(header1, bedlines1, ignore_attr)
        bed_str2 = BedLineParser.to_str(header2, bedlines2, ignore_attr)

    # Find replace the prefix
    if ignore_primer_prefix:
        bl1_pn_prefixs = {bl.amplicon_prefix for bl in bedlines1}
        for bl1_pn_prefix in bl1_pn_prefixs:
            bed_str1 = bed_str1.replace(bl1_pn_prefix, PLACEHOLDER_STR)

        bl2_pn_prefixs = {bl.amplicon_prefix for bl in bedlines2}
        for bl2_pn_prefix in bl2_pn_prefixs:
            bed_str2 = bed_str2.replace(bl2_pn_prefix, PLACEHOLDER_STR)

    return (bed_str1, bed_str2)


def diff_primernames(
    bedlines1: list[BedLine], bedlines2: list[BedLine]
) -> tuple[set[str], set[str]]:
    """Returns difference set operations on the primernames.

    Calculates the set difference of primer names between two lists of BedLines.

    Args:
        bedlines1 (list[BedLine]): The first list of BedLine objects.
        bedlines2 (list[BedLine]): The second list of BedLine objects.

    Returns:
        tuple[set[str], set[str]]: A tuple containing two sets:
            - Primer names present in bedlines1 but not in bedlines2.
            - Primer names present in bedlines2 but not in bedlines1.
    """
    # Get primernames present in each.
    primernames1 = set([bl.primername for bl in bedlines1])
    primernames2 = set([bl.primername for bl in bedlines2])

    # Returns primernames present in X but not Y.
    return primernames1.difference(primernames2), primernames2.difference(primernames1)


def diff_sequence(
    bedlines1: list[BedLine], bedlines2: list[BedLine]
) -> tuple[set[str], set[str]]:
    """Returns difference set operations on the sequences.

    Calculates the set difference of sequences between two lists of BedLines.

    Args:
        bedlines1 (list[BedLine]): The first list of BedLine objects.
        bedlines2 (list[BedLine]): The second list of BedLine objects.

    Returns:
        tuple[set[str], set[str]]: A tuple containing two sets:
            - Sequences present in bedlines1 but not in bedlines2.
            - Sequences present in bedlines2 but not in bedlines1.
    """
    # Get primernames present in each.
    sequences1 = set([bl.sequence for bl in bedlines1])
    sequences2 = set([bl.sequence for bl in bedlines2])

    # Returns primernames present in X but not Y.
    return sequences1.difference(sequences2), sequences2.difference(sequences1)


def ndiff_bedlines(
    bedlines1: list[BedLine],
    bedlines2: list[BedLine],
    header1: Optional[list[str]] = None,
    header2: Optional[list[str]] = None,
    ignore_order: bool = True,
    ignore_attr: bool = False,
    ignore_header: bool = False,
    ignore_no_diff: bool = False,
    ignore_primer_prefix: bool = False,
):
    """Generates a difference report between two sets of bedlines using difflib.ndiff.

    Compares two lists of BedLines and produces a generator of differences, similar to
    the `ndiff` tool. Can optionally ignore order, attributes, and headers during comparison.

    Args:
        bedlines1 (list[BedLine]): The first list of BedLine objects.
        bedlines2 (list[BedLine]): The second list of BedLine objects.
        header1 (Optional[list[str]] , optional): Headers for the first bedfile. Defaults to None.
        header2 (Optional[list[str]] , optional): Headers for the second bedfile. Defaults to None.
        ignore_order (bool, optional): If True, sorts bedlines before comparison. Defaults to True.
        ignore_attr (bool, optional): If True, excludes attributes from comparison. Defaults to False.
        ignore_header (bool, optional): If True, excludes headers from comparison. Defaults to False.
        ignore_no_diff (bool, optional): If True, filters out lines that are identical (starting with "  "). Defaults to False.
        ignore_primer_prefix (bool, optional): If True, excludes the primername prefix from the string representation. Defaults to False.

    Returns:
        Iterator[str]: A generator yielding difference lines (strings).
    """
    nbed_str1, nbed_str2 = create_normalised_bedfile_str(
        bedlines1=bedlines1,
        bedlines2=bedlines2,
        header1=header1,
        header2=header2,
        ignore_attr=ignore_attr,
        ignore_order=ignore_order,
        ignore_header=ignore_header,
        ignore_primer_prefix=ignore_primer_prefix,
    )

    diff_generator = ndiff(
        nbed_str1.splitlines(keepends=True),
        nbed_str2.splitlines(keepends=True),
    )
    if ignore_no_diff:
        return (line for line in diff_generator if not line.startswith("  "))
    return diff_generator


def unified_diff_bedlines(
    bedlines1: list[BedLine],
    bedlines2: list[BedLine],
    header1: Optional[list[str]] = None,
    header2: Optional[list[str]] = None,
    ignore_order: bool = True,
    ignore_attr: bool = False,
    ignore_header: bool = False,
    ignore_primer_prefix: bool = False,
):
    """Generates a unified difference report between two sets of bedlines.

    Compares two lists of BedLines and produces a generator of differences in unified diff format.
    Can optionally ignore order, attributes, and headers during comparison.

    Args:
        bedlines1 (list[BedLine]): The first list of BedLine objects.
        bedlines2 (list[BedLine]): The second list of BedLine objects.
        header1 (Optional[list[str]] , optional): Headers for the first bedfile. Defaults to None.
        header2 (Optional[list[str]] , optional): Headers for the second bedfile. Defaults to None.
        ignore_order (bool, optional): If True, sorts bedlines before comparison. Defaults to True.
        ignore_attr (bool, optional): If True, excludes attributes from comparison. Defaults to False.
        ignore_header (bool, optional): If True, excludes headers from comparison. Defaults to False.
        ignore_primer_prefix (bool, optional): If True, excludes the primername prefix from the string representation. Defaults to False.

    Returns:
        Iterator[str]: A generator yielding unified diff lines (strings).
    """
    # Try and match on pn.

    nbed_str1, nbed_str2 = create_normalised_bedfile_str(
        bedlines1=bedlines1,
        bedlines2=bedlines2,
        header1=header1,
        header2=header2,
        ignore_attr=ignore_attr,
        ignore_order=ignore_order,
        ignore_header=ignore_header,
        ignore_primer_prefix=ignore_primer_prefix,
    )

    return unified_diff(
        nbed_str1.splitlines(keepends=True),
        nbed_str2.splitlines(keepends=True),
        n=0,
        fromfile="bedlines1",
        tofile="bedlines2",
    )
