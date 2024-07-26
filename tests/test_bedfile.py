import pathlib
import unittest

from primalbedtools.bedfiles import (
    BedLine,
    create_bedfile_str,
    create_bedline,
    group_by_amplicon_number,
    group_by_chrom,
    group_by_strand,
    read_bedfile,
    write_bedfile,
)

TEST_BEDFILE = pathlib.Path(__file__).parent / "test.bed"


class TestBedLine(unittest.TestCase):
    def test_bedline_create(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 1)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(bedline.ipool, 0)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n",
        )


class TestCreateBedline(unittest.TestCase):
    def test_create_bedline(self):
        bedline = create_bedline(
            ["chr1", "100", "200", "scheme_1_LEFT", "1", "+", "ACGT"]
        )
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")


class TestReadBedfile(unittest.TestCase):
    def test_read_bedfile(self):
        headers, bedlines = read_bedfile(TEST_BEDFILE)
        self.assertEqual(
            headers, ["# artic-bed-version v3.0", "# artic-sars-cov-2 / 400 / v5.3.2"]
        )

        self.assertEqual(len(bedlines), 6)
        self.assertEqual(bedlines[0].chrom, "MN908947.3")


class TestCreateBedfileStr(unittest.TestCase):
    bedline = BedLine(
        chrom="chr1",
        start=100,
        end=200,
        primername="scheme_1_LEFT",
        pool=1,
        strand="+",
        sequence="ACGT",
    )

    def test_create_bedfile_str(self):
        bedfile_str = create_bedfile_str(["#header1"], [self.bedline])
        self.assertEqual(
            bedfile_str, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )

    def test_create_bedfile_str_no_header(self):
        bedfile_str = create_bedfile_str([], [self.bedline])
        self.assertEqual(bedfile_str, "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n")

    def test_create_bedfile_str_malformed_header(self):
        bedfile_str = create_bedfile_str(["header1"], [self.bedline])
        self.assertEqual(
            bedfile_str, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )


class TestWriteBedfile(unittest.TestCase):
    output_bed_path = pathlib.Path(__file__).parent / "test_output.bed"

    def test_write_bedfile(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        write_bedfile(self.output_bed_path, ["#header1"], [bedline])
        with open("test_output.bed") as f:
            content = f.read()
        self.assertEqual(
            content, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )

    def tearDown(self) -> None:
        self.output_bed_path.unlink()
        super().tearDown()


class TestGroupByChrom(unittest.TestCase):
    def test_group_by_chrom(self):
        bedline1 = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        bedline2 = BedLine(
            chrom="chr2",
            start=150,
            end=250,
            primername="scheme_2_LEFT",
            pool=2,
            strand="+",
            sequence="ACGT",
        )
        grouped = group_by_chrom([bedline1, bedline2])
        self.assertEqual(len(grouped), 2)
        self.assertEqual(len(grouped["chr1"]), 1)
        self.assertEqual(len(grouped["chr2"]), 1)
        self.assertEqual(grouped["chr1"][0].chrom, "chr1")
        self.assertEqual(grouped["chr2"][0].chrom, "chr2")

    def test_group_by_chrom_empty(self):
        grouped = group_by_chrom([])
        self.assertEqual(len(grouped), 0)

    def test_group_by_chrom_file(self):
        headers, bedlines = read_bedfile(TEST_BEDFILE)
        grouped = group_by_chrom(bedlines)
        self.assertEqual(len(grouped), 1)
        self.assertEqual(len(grouped["MN908947.3"]), 6)


class TestGroupByAmpliconNumber(unittest.TestCase):
    def test_group_by_amplicon_number(self):
        bedline1 = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        bedline2 = BedLine(
            chrom="chr1",
            start=150,
            end=250,
            primername="scheme_1_LEFT",
            pool=2,
            strand="+",
            sequence="ACGT",
        )
        grouped = group_by_amplicon_number([bedline1, bedline2])
        self.assertEqual(len(grouped), 1)
        self.assertEqual(len(grouped[1]), 2)
        self.assertEqual(grouped[1][0].chrom, "chr1")
        self.assertEqual(grouped[1][1].chrom, "chr1")

    def test_group_by_amplicon_number_file(self):
        headers, bedlines = read_bedfile(TEST_BEDFILE)
        grouped = group_by_amplicon_number(bedlines)
        self.assertEqual(len(grouped), 3)

        # Check for correct ampliconnumber
        for amplicon_number, bedlines in grouped.items():
            for bedline in bedlines:
                self.assertEqual(bedline.amplicon_number, amplicon_number)

        # Check for correct primernames
        self.assertEqual(
            {bl.primername for bl in grouped[1]},
            {"SARS-CoV-2_1_LEFT_1", "SARS-CoV-2_1_RIGHT_1"},
        )
        self.assertEqual(
            {bl.primername for bl in grouped[2]},
            {"SARS-CoV-2_2_LEFT_0", "SARS-CoV-2_2_RIGHT_0"},
        )
        self.assertEqual(
            {bl.primername for bl in grouped[3]},
            {"SARS-CoV-2_3_LEFT_1", "SARS-CoV-2_3_RIGHT_0"},
        )


class TestGroupByStrand(unittest.TestCase):
    def test_group_by_strand(self):
        bedline1 = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        bedline2 = BedLine(
            chrom="chr1",
            start=150,
            end=250,
            primername="scheme_1_RIGHT",
            pool=2,
            strand="-",
            sequence="ACGT",
        )
        grouped = group_by_strand([bedline1, bedline2])
        self.assertEqual(len(grouped), 2)
        self.assertEqual(len(grouped["+"]), 1)
        self.assertEqual(len(grouped["-"]), 1)
        self.assertEqual(grouped["+"], [bedline1])
        self.assertEqual(grouped["-"], [bedline2])

    def test_group_by_strand_file(self):
        headers, bedlines = read_bedfile(TEST_BEDFILE)
        grouped = group_by_strand(bedlines)
        self.assertEqual(len(grouped), 2)
        self.assertEqual(len(grouped["+"]), 3)
        self.assertEqual(len(grouped["-"]), 3)


if __name__ == "__main__":
    unittest.main()
