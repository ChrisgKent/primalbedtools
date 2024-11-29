import io
import unittest

from primalbedtools.bedfiles import BedLine
from primalbedtools.fasta import read_fasta
from primalbedtools.remap import (
    create_mapping_list,
    remap,
)


class TestMappingList(unittest.TestCase):
    def test_mapping_list_simple(self):
        """
        Test mapping list creation with no gaps
        """
        fbedline = BedLine("chr1", 10, 20, "test_1_LEFT_1", 1, "+", "ATCGATCGAA")
        rbedline = BedLine(
            "chr1", 30, 40, "test_1_RIGHT_1", 1, "-", "ATCGATCGAA"
        )  # seq normally rc

        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNATCGATCGAA\n"
            ">chr2\n"
            "----------ATCGATCGAANNNNNNNNNNATCGATCGAA\n"
        )
        msa = read_fasta(fasta_io)
        msa_to_genome, from_index_to_msa_index = create_mapping_list(
            msa, "chr1", "chr2"
        )

        # Check the reference to msa mapping
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[fbedline.start]], fbedline.start
        )
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[fbedline.end]], fbedline.end
        )
        # Check the mapping to the new reference
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[fbedline.start]], 0)
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[fbedline.end]], 10)

        # Test the ref sequence is a match
        ref_seq_str = msa["chr1"].replace("-", "")
        self.assertEqual(
            ref_seq_str[
                msa_to_genome[0][
                    from_index_to_msa_index[fbedline.start]
                ] : msa_to_genome[0][from_index_to_msa_index[fbedline.end]]
            ].replace("-", ""),
            fbedline.sequence,
        )

        # Check for the reverse primer
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[rbedline.start]], rbedline.start
        )
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[rbedline.end]], rbedline.end
        )
        # Check the mapping to the new reference
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[rbedline.start]], 20)
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[rbedline.end]], 30)
        # Test the ref sequence is a match
        new_seq_str = msa["chr2"].replace("-", "")
        self.assertEqual(
            new_seq_str[
                msa_to_genome[1][
                    from_index_to_msa_index[rbedline.start]
                ] : msa_to_genome[1][from_index_to_msa_index[rbedline.end]]
            ].replace("-", ""),
            rbedline.sequence,
        )

    def test_mapping_list_no_id(self):
        """
        Test mapping list creation when ID is not found
        """
        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNTTCGATCGATN\n"
            ">chr2\n"
            "----------ATCGATCGAANNNNNNNNNNTTCGATCGATN\n"
        )
        msa = read_fasta(fasta_io)
        with self.assertRaises(ValueError):
            create_mapping_list(msa, "chr1", "chr3")

    def test_mapping_list_dif_len(self):
        """
        Test mapping list creation when the lengths are different
        """
        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNTTCGATCGATN\n"
            ">chr2\n"
            "----------ATCGATCGAANNNNNNNNNNTTCGATCGATNN\n"
        )
        msa = read_fasta(fasta_io)
        with self.assertRaises(ValueError):
            create_mapping_list(msa, "chr1", "chr2")

    def test_mapping_list_gaps(self):
        """
        Tests how mapping list handles gaps
        """
        fbedline = BedLine("chr1", 10, 20, "test_1_LEFT_1", 1, "+", "ATCGATCGAA")
        rbedline = BedLine("chr1", 30, 40, "test_1_RIGHT_1", 1, "-", "ATCGATCGAA")

        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCG-ATC--GAANNNNNNNNNNATCGATCGAA\n"
            ">chr2\n"
            "----------ATCG-ATC--GAANNNNNNNNNNATCGATCGAA\n"
        )
        msa = read_fasta(fasta_io)
        msa_to_genome, from_index_to_msa_index = create_mapping_list(
            msa, "chr1", "chr2"
        )

        # Check the reference to msa mapping
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[fbedline.start]], fbedline.start
        )
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[fbedline.end]], fbedline.end
        )
        # Check the mapping to the new reference
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[fbedline.start]], 0)
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[fbedline.end]], 10)

        # Test the ref sequence is a match
        ref_seq_str = msa["chr1"].replace("-", "")
        self.assertEqual(
            ref_seq_str[
                msa_to_genome[0][
                    from_index_to_msa_index[fbedline.start]
                ] : msa_to_genome[0][from_index_to_msa_index[fbedline.end]]
            ].replace("-", ""),
            fbedline.sequence,
        )

        # Check for the reverse primer
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[rbedline.start]], rbedline.start
        )
        self.assertEqual(
            msa_to_genome[0][from_index_to_msa_index[rbedline.end]], rbedline.end
        )
        # Check the mapping to the new reference
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[rbedline.start]], 20)
        self.assertEqual(msa_to_genome[1][from_index_to_msa_index[rbedline.end]], 30)
        # Test the ref sequence is a match
        new_seq_str = msa["chr2"].replace("-", "")
        self.assertEqual(
            new_seq_str[
                msa_to_genome[1][
                    from_index_to_msa_index[rbedline.start]
                ] : msa_to_genome[1][from_index_to_msa_index[rbedline.end]]
            ].replace("-", ""),
            rbedline.sequence,
        )


class Testremap(unittest.TestCase):
    def test_remap_simple(self):
        """
        Test perfect remapping with no gaps
        """
        # Bedline to remap
        fbedline = BedLine("chr1", 10, 20, "test_1_LEFT_1", 1, "+", "ATCGATCGAA")
        rbedline = BedLine("chr1", 30, 40, "test_1_RIGHT_1", 1, "-", "TTCGATCGAT")

        # MSAs to remap to
        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNTTCGATCGATN\n"
            ">chr2\n"
            "----------ATCGATCGAANNNNNNNNNNTTCGATCGATN\n"
        )
        msa = read_fasta(fasta_io)
        # Check primer maps corrected
        self.assertEqual(
            (msa[fbedline.chrom][fbedline.start : fbedline.end]),
            fbedline.sequence,
        )

        # Remap
        remap("chr1", "chr2", [fbedline, rbedline], msa)

        # Check primer maps corrected
        self.assertEqual(fbedline.start, 0)
        self.assertEqual(fbedline.end, 10)
        self.assertEqual(rbedline.start, 20)
        self.assertEqual(rbedline.end, 30)

        # Check chromname is updated
        self.assertEqual(fbedline.chrom, "chr2")
        self.assertEqual(rbedline.chrom, "chr2")

    def test_remap_full_gap(self):
        """
        Test remapping when primer is not present
        """
        # Bedline to remap
        fbedline = BedLine("chr1", 10, 20, "test_1_LEFT_1", 1, "+", "ATCGATCGAA")
        rbedline = BedLine("chr1", 30, 40, "test_1_RIGHT_1", 1, "-", "ATCGATCGAA")

        # MSAs to remap to
        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNTTCGATCGATN\n"
            ">chr2\n"
            "--------------------NNNNNNNNNNN----------\n"
        )
        msa = read_fasta(fasta_io)

        # Remaps
        remap("chr1", "chr2", [fbedline], msa)

        # Check nothing changes
        self.assertEqual(fbedline.start, 10)
        self.assertEqual(fbedline.end, 20)
        self.assertEqual(fbedline.chrom, "chr1")

        self.assertEqual(rbedline.start, 30)
        self.assertEqual(rbedline.end, 40)
        self.assertEqual(rbedline.chrom, "chr1")

    def test_remap_no_3p_gap_at_edge(self):
        """
        Test remapping when primer doesn't have a 3' gap. At the limits of the MSA
        """
        # Bedline to remap
        fbedline = BedLine("chr1", 10, 20, "test_1_LEFT_1", 1, "+", "ATCGATCGAA")
        rbedline = BedLine(
            "chr1", 30, 40, "test_1_RIGHT_1", 1, "-", "ATCGATCGAA"
        )  # seq normally rc

        # MSAs to remap to
        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNATCGATCGAA\n"
            ">chr2\n"
            "----------ATC-ATCGAANNNNNNNNNNATCGATCGAA\n"
        )
        msa = read_fasta(fasta_io)

        # Remaps
        remap("chr1", "chr2", [fbedline, rbedline], msa)

        print(fbedline.to_bed())

        # Check primer maps corrected
        self.assertEqual(fbedline.start, 0)
        self.assertEqual(fbedline.end, 9)

        self.assertEqual(rbedline.start, 19)
        self.assertEqual(rbedline.end, 29)

        # Check chromname is updated
        self.assertEqual(fbedline.chrom, "chr2")
        self.assertEqual(rbedline.chrom, "chr2")

    def test_remap_with_3p_indel(self):
        """
        Test remapping when primer has a 3' gap
        """
        # Bedline to remap
        fbedline = BedLine("chr1", 10, 20, "test_1_LEFT_1", 1, "+", "ATCGATCGAA")
        rbedline = BedLine("chr1", 30, 40, "test_1_RIGHT_1", 1, "-", "ATCGATCGAA")

        # MSAs to remap to
        fasta_io = io.StringIO(
            ">chr1\n"
            "NNNNNNNNNNATCGATCGAANNNNNNNNNNATCGATCGAAATCGATCGATCGATCGATCG\n"
            ">chr2\n"
            "NNNNNNNNNNATCGATCGA-NNNNNNNNNN-TCGATCGAAATCGATCGATCGATCGATCG\n"
        )
        msa = read_fasta(fasta_io)

        # Remaps
        remap("chr1", "chr2", [fbedline, rbedline], msa)

        # Check primer maps corrected
        self.assertEqual(fbedline.start, 10)
        self.assertEqual(fbedline.end, 20)

        self.assertEqual(rbedline.start, 28)
        self.assertEqual(rbedline.end, 38)


if __name__ == "__main__":
    unittest.main()
