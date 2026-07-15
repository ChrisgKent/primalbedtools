# To keep deps low here is a simple fasta parser

from io import TextIOBase
from typing import Union


def read_fasta(fasta_file: Union[str, TextIOBase]) -> dict[str, str]:
    """
    Read a fasta file and return a dictionary with the sequence name as the key and the sequence as the value.
    """
    sequences = {}

    if isinstance(fasta_file, str):
        handle = open(fasta_file)
    else:
        handle = fasta_file

    seq_name = None
    with handle as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seq_name = line[1:].split()[0]
                if seq_name in sequences:
                    raise ValueError(
                        f"Duplicate sequence name: {seq_name} (line {lineno})"
                    )
                sequences[seq_name] = []
            else:
                if seq_name is None:
                    raise ValueError(
                        f"Sequence data found before any header (line {lineno}): {line[:50]!r}"
                    )
                sequences[seq_name].append(line)  # type: ignore

    # Avoid str concatenation
    for seq_name, seq in sequences.items():
        sequences[seq_name] = "".join(seq)

    return sequences
