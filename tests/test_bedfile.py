import copy
import pathlib
import random
import unittest

from primalbedtools.bedfiles import (
    BedLine,
    BedLineParser,
    PrimerNameVersion,
    create_bedfile_str,
    create_bedline,
    downgrade_primernames,
    group_by_amplicon_number,
    group_by_chrom,
    group_by_strand,
    read_bedfile,
    sort_bedlines,
    update_primernames,
    version_primername,
    write_bedfile,
)

TEST_BEDFILE = pathlib.Path(__file__).parent / "test.bed"
TEST_V2_BEDFILE = pathlib.Path(__file__).parent / "test.v2.bed"


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

    def test_invalid_bedline(self):
        # Fake primername should raise ValueError
        with self.assertRaises(ValueError):
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="fake_primername",
                pool=1,
                strand="+",
                sequence="ACGT",
            )
        # 0-based pool should raise ValueError
        with self.assertRaises(ValueError):
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="scheme_1_LEFT",
                pool=0,
                strand="+",
                sequence="ACGT",
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

    def test_read_v2_bedfile(self):
        headers, bedlines = read_bedfile(TEST_V2_BEDFILE)
        # Check for empty headers
        self.assertEqual(headers, [])

        # Check correct number of bedlines and chrom
        self.assertEqual(len(bedlines), 10)
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
        with open(self.output_bed_path) as f:
            content = f.read()
        self.assertEqual(
            content, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )

    def tearDown(self) -> None:
        self.output_bed_path.unlink(missing_ok=True)
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
        """
        Tests grouping by chrom for a v3 (default) bedfile
        """
        headers, bedlines = read_bedfile(TEST_BEDFILE)
        grouped = group_by_chrom(bedlines)
        self.assertEqual(len(grouped), 1)
        self.assertEqual(len(grouped["MN908947.3"]), 6)

    def test_group_by_chrom_v2_file(self):
        """
        Tests grouping by chrom for a v2 bedfile
        """
        headers, bedlines = read_bedfile(TEST_V2_BEDFILE)
        grouped = group_by_chrom(bedlines)
        self.assertEqual(len(grouped), 1)
        self.assertEqual(len(grouped["MN908947.3"]), 10)


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


class TestBedLineParser(unittest.TestCase):
    OUTFILE = pathlib.Path(__file__).parent / "test_bedline_parser_to_file.bed"

    def test_bedline_parser_from_file(self):
        headers, bedlines = BedLineParser.from_file(TEST_BEDFILE)
        self.assertEqual(
            headers, ["# artic-bed-version v3.0", "# artic-sars-cov-2 / 400 / v5.3.2"]
        )

        self.assertEqual(len(bedlines), 6)
        self.assertEqual(bedlines[0].chrom, "MN908947.3")

    def test_bedline_parser_from_str(self):
        with open(TEST_BEDFILE) as f:
            bedfile_str = f.read()
        headers, bedlines = BedLineParser.from_str(bedfile_str)
        self.assertEqual(
            headers, ["# artic-bed-version v3.0", "# artic-sars-cov-2 / 400 / v5.3.2"]
        )

        self.assertEqual(len(bedlines), 6)
        self.assertEqual(bedlines[0].chrom, "MN908947.3")

    def test_bedline_parser_to_str(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        bedfile_str = BedLineParser.to_str(["#header1"], [bedline])
        self.assertEqual(
            bedfile_str, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )

    def test_bedline_parser_to_file(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        BedLineParser.to_file(self.OUTFILE, ["#header1"], [bedline])
        with open(self.OUTFILE) as f:
            content = f.read()
        self.assertEqual(
            content, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )

    def tearDown(self) -> None:
        self.OUTFILE.unlink(missing_ok=True)
        super().tearDown()


class TestModifyBedLines(unittest.TestCase):
    v2_bedlines = BedLineParser.from_file(TEST_V2_BEDFILE)[1]

    def test_update_primername_simple(self):
        local_v2_bedlines = copy.deepcopy(self.v2_bedlines)
        old_bednames = {bl.primername for bl in local_v2_bedlines}

        # Update primernames
        v3_bedlines = update_primernames(local_v2_bedlines)
        new_bednames = {bl.primername for bl in v3_bedlines}

        # Check same length
        self.assertEqual(len(old_bednames), len(new_bednames))

        # Check for correct primernames
        self.assertEqual(
            new_bednames,
            {
                "SARS-CoV-2_1_LEFT_1",
                "SARS-CoV-2_1_RIGHT_1",
                "SARS-CoV-2_2_LEFT_1",
                "SARS-CoV-2_2_RIGHT_1",
                "SARS-CoV-2_3_LEFT_1",
                "SARS-CoV-2_3_RIGHT_1",
                "SARS-CoV-2_4_LEFT_1",
                "SARS-CoV-2_4_RIGHT_1",
                "SARS-CoV-2_5_LEFT_1",
                "SARS-CoV-2_5_RIGHT_1",
            },
        )

        # Check all the new names are v2
        self.assertTrue(
            all(
                version_primername(name) == PrimerNameVersion.V2
                for name in new_bednames
            )
        )

    def test_update_primername_alt(self):
        bedlines = [
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="test_1_LEFT_alt",
                pool=1,
                strand="+",
                sequence="ACGT",
            ),
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="test_1_LEFT",
                pool=1,
                strand="+",
                sequence="ACGT",
            ),
        ]
        new_bedlines = update_primernames(bedlines)
        new_primername = {bl.primername for bl in new_bedlines}
        self.assertEqual(new_primername, {"test_1_LEFT_1", "test_1_LEFT_2"})

    def test_downgrade_primername(self):
        bedlines = [
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="test_1_LEFT",
                pool=1,
                strand="+",
                sequence="ACGT",
            ),
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="test_1_LEFT_alt",
                pool=1,
                strand="+",
                sequence="ACGT",
            ),
        ]
        new_bedlines = downgrade_primernames(bedlines)
        new_primername = {bl.primername for bl in new_bedlines}
        self.assertEqual(new_primername, {"test_1_LEFT", "test_1_LEFT_alt1"})

    def test_sort_bedlines(self):
        # Read in a bedfile
        headers, bedlines = BedLineParser.from_file(TEST_BEDFILE)

        # Randomly shuffle the bedlines
        random.seed(100)
        random_bedlines = random.sample(bedlines, len(bedlines))

        # Sort the bedlines
        sorted_bedlines = sort_bedlines(random_bedlines)

        # Check that the bedlines are sorted
        self.assertEqual(sorted_bedlines, bedlines)


if __name__ == "__main__":
    unittest.main()
