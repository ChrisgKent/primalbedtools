AMBIGUOUS_DNA_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    "-": "-",
}

ALL_DNA_WITH_N: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "ACGT",
}


def expand_ambiguous_bases(seq: str) -> list[str]:
    """
    Expand ambiguous bases in a DNA sequence to all possible combinations.
    """
    if not seq:
        return []
    # Initialize the list with the first character
    expanded = [b for b in ALL_DNA_WITH_N[seq[0]]]
    # Iterate through the sequence and expand ambiguous bases
    for base in seq[1:]:
        new_expanded = []
        for current in expanded:
            for option in ALL_DNA_WITH_N[base]:  # raise KeyError
                new_expanded.append(
                    current + option
                )  # slow but not performance critical

        expanded = new_expanded
    return expanded


def rc_seq(seq: str) -> str:
    """
    Reverse complement a DNA sequence.
    """
    return complement_seq(seq[::-1])


def complement_seq(seq: str) -> str:
    """
    Complement a DNA sequence.
    """
    return "".join(AMBIGUOUS_DNA_COMPLEMENT[base] for base in seq)


def strip_all_white_space(s: str) -> str:
    """Strips any whitespace"""
    return "".join(s.split())
