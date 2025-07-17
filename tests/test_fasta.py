import pathlib
import unittest
from io import StringIO

from primalbedtools.fasta import read_fasta

FASTA_PATH = pathlib.Path(__file__).parent / "inputs/msa.input.fasta"


class TestFasta(unittest.TestCase):
    def test_read_fasta_from_handle(self):
        fasta_io = StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCG-ATC--GAANNNNNNNNNNATCGATCGAA\n"
            ">chr2\n"
            "----------ATCG-ATC--GAANNNNNNNNNNATCGATCGAA\n"
        )
        msa = read_fasta(fasta_io)
        self.assertEqual(msa["chr1"], "NNNNNNNNNNATCG-ATC--GAANNNNNNNNNNATCGATCGAA")
        self.assertEqual(msa["chr2"], "----------ATCG-ATC--GAANNNNNNNNNNATCGATCGAA")

    def test_read_fasta_from_file(self):
        fasta_path = FASTA_PATH.resolve()
        msa = read_fasta(str(fasta_path))

        self.assertEqual(msa.keys(), {"seq1", "seq2"})

        self.assertEqual(
            msa["seq1"],
            "ATCGATCGATCATCGATCGATCGTAGCTAGCAYCGCTAGCTAGCGATCGATCGCAYTGCACCCAACCATGTACCGTCGAGTTA",
        )

        self.assertEqual(msa["seq2"], "ATCGATCGATCATCGATCGAT")
