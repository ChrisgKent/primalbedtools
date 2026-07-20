import unittest
from io import StringIO

from primalbedtools.fasta import read_fasta
from tests.infiles import FASTA_PATH


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

    def test_leading_blank_lines_are_ignored(self):
        fasta_io = StringIO("\n\n>seq1\nACGT\n")
        msa = read_fasta(fasta_io)
        self.assertEqual(msa["seq1"], "ACGT")

    def test_sequence_before_header_raises(self):
        fasta_io = StringIO("ACGT\n>seq1\nACGT\n")
        # The message should name the offending line number and content.
        with self.assertRaisesRegex(ValueError, r"before any header \(line 1\)"):
            read_fasta(fasta_io)

    def test_duplicate_sequence_name_raises_with_line(self):
        fasta_io = StringIO(">seq1\nACGT\n>seq1\nTTTT\n")
        # The duplicate header is on line 3 and should be named.
        with self.assertRaisesRegex(
            ValueError, r"Duplicate sequence name: seq1 \(line 3\)"
        ):
            read_fasta(fasta_io)
