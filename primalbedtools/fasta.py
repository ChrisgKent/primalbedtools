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

    with handle as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line[1:].split()[0]
                if seq_name in sequences:
                    raise ValueError(f"Duplicate sequence name: {seq_name}")
                sequences[seq_name] = []
            else:
                sequences[seq_name].append(line)

    # Avoid str concatenation
    for seq_name, seq in sequences.items():
        sequences[seq_name] = "".join(seq)

    return sequences
