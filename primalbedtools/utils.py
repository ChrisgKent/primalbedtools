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
