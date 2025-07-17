import pathlib
import unittest

from primalbedtools.amplicons import Amplicon
from primalbedtools.bedfiles import BedLine
from primalbedtools.validate import (
    do_pp_ol,
    validate,
    validate_primerbed,
    validate_ref_and_bed,
)

FASTA_PATH = pathlib.Path(__file__).parent / "inputs/msa.input.fasta"


class TestValidate(unittest.TestCase):
    def setUp(self) -> None:
        self.fbl1 = BedLine(
            chrom="chr1",
            start=100,
            end=120,
            primername="scheme_1_LEFT_0",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        self.rbl1 = BedLine(
            chrom="chr1",
            start=200,
            end=220,
            primername="scheme_1_RIGHT_0",
            pool=1,
            strand="-",
            sequence="ACGT",
        )
        self.fbl2 = BedLine(
            chrom="chr1",
            start=150,
            end=160,
            primername="scheme_2_LEFT_0",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        self.rbl2 = BedLine(
            chrom="chr1",
            start=250,
            end=260,
            primername="scheme_2_RIGHT_0",
            pool=1,
            strand="-",
            sequence="ACGT",
        )

        self.fasta_path = str(FASTA_PATH.absolute())
        return super().setUp()

    def test_do_pp_ol(self):
        pp1 = Amplicon([self.fbl1], [self.rbl1])
        pp2 = Amplicon([self.fbl2], [self.rbl2])

        # Detect ol
        self.assertTrue(do_pp_ol(pp1, pp2))

        # Change region
        self.fbl2.start = 1000
        self.fbl2.end = 1010

        self.rbl2.start = 1200
        self.rbl2.start = 1210
        # Check no ol
        self.assertFalse(do_pp_ol(pp1, pp2))

    def test_validate_primerbed_ol(self):
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # overlap
        with self.assertRaises(ValueError) as _:
            validate_primerbed(bls)

        # Change region
        self.fbl2.start = 1000
        self.fbl2.end = 1010

        self.rbl2.start = 1200
        self.rbl2.start = 1210

        # overlap
        validate_primerbed(bls)

    def test_validate_primerbed_diff_pool(self):
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # pool
        with self.assertRaises(ValueError) as _:
            validate_primerbed(bls)

        # Change pool
        self.fbl2.pool = 2
        self.rbl2.pool = 2

        # overlap
        validate_primerbed(bls)

    def test_validate_primerbed_diff_chrom(self):
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # Change chrom
        self.fbl2.chrom = "test"
        self.rbl2.chrom = "test"

        # overlap
        validate_primerbed(bls)

    def test_validate_primerbed_no_left(self):
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # Change pool
        self.fbl2.pool = 2
        self.rbl2.pool = 2

        # check pass
        validate_primerbed(bls)

        bls = bls[1:]

        # check fail
        with self.assertRaises(ValueError) as _:
            validate_primerbed(bls)

    def test_validate_primerbed_no_right(self):
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # Change pool
        self.fbl2.pool = 2
        self.rbl2.pool = 2

        # check pass
        validate_primerbed(bls)

        bls = bls[:-1]

        # check fail
        with self.assertRaises(ValueError) as _:
            validate_primerbed(bls)

    def test_validate_ref_and_bed(self):
        """
        Tests matching chrom passes
        """
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # Set chromnames
        self.fbl1.chrom = "seq1"
        self.rbl1.chrom = "seq1"

        self.fbl2.chrom = "seq2"
        self.rbl2.chrom = "seq2"

        # Check valid
        validate_ref_and_bed(bls, self.fasta_path)

    def test_validate_ref_and_bed_extra_fasta(self):
        """
        Chrom in fasta not in bedfile
        """
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # Set chromnames
        self.fbl1.chrom = "seq1"
        self.rbl1.chrom = "seq1"
        self.fbl2.chrom = "seq1"
        self.rbl2.chrom = "seq1"

        with self.assertRaises(ValueError) as cm:
            validate_ref_and_bed(bls, self.fasta_path)

        # Check correct chrom comes up
        self.assertIn("seq2", str(cm.exception))
        self.assertIn(
            "chroms in reference.fasta are not in primer.bed", str(cm.exception)
        )

    def test_validate_ref_and_bed_extra_bed(self):
        """
        Chrom in fasta not in bedfile
        """
        bls = [self.fbl1, self.rbl1, self.fbl2, self.rbl2]

        # Set chromnames
        self.fbl1.chrom = "seq10"
        self.rbl1.chrom = "seq10"
        self.fbl2.chrom = "seq10"
        self.rbl2.chrom = "seq10"

        with self.assertRaises(ValueError) as cm:
            validate_ref_and_bed(bls, self.fasta_path)

        # Check correct chrom comes up
        self.assertIn(self.fbl1.chrom, str(cm.exception))
        self.assertIn(
            "chroms in primer.bed are not in reference.fasta", str(cm.exception)
        )

    def test_validate(self):
        """Test validate with files"""
        validate("tests/inputs/primer.bed", "tests/inputs/reference.fasta")


if __name__ == "__main__":
    unittest.main()
