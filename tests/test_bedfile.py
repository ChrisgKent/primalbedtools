import copy
import pathlib
import random
import unittest

from primalbedtools.bedfiles import (
    PRIMER_WEIGHT_KEY,
    BedLine,
    BedLineParser,
    PrimerClass,
    PrimerNameVersion,
    Strand,
    create_bedfile_str,
    create_bedline,
    create_primer_attributes_str,
    downgrade_primernames,
    expand_bedlines,
    group_by_amplicon_number,
    group_by_chrom,
    group_by_class,
    group_by_pool,
    group_by_strand,
    lr_string_to_strand_char,
    merge_primers,
    parse_headers_to_dict,
    parse_primer_attributes_str,
    primer_class_str_to_enum,
    read_bedfile,
    sort_bedlines,
    strand_char_to_primer_class_str,
    update_primernames,
    version_primername,
    write_bedfile,
)

TEST_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.bed"
TEST_V2_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.v2.bed"
TEST_WEIGHTS_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.weights.bed"
TEST_WEIGHTS_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.weights.bed"
TEST_ATTRIBUTES_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.attributes.bed"
TEST_PROBE_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.probe.bed"


class TestValidationFuncs(unittest.TestCase):
    def test_string_to_strand_char(self):
        # Check expected
        self.assertEqual(lr_string_to_strand_char("LEFT"), "+")
        self.assertEqual(lr_string_to_strand_char("RIGHT"), "-")

        # Check unexpected
        with self.assertRaises(ValueError):
            lr_string_to_strand_char("")

    def test_primer_class_str_to_enum(self):
        self.assertEqual(primer_class_str_to_enum("LEFT"), PrimerClass.LEFT)
        self.assertEqual(primer_class_str_to_enum("RIGHT"), PrimerClass.RIGHT)
        self.assertEqual(primer_class_str_to_enum("PROBE"), PrimerClass.PROBE)

        # Check unexpected
        with self.assertRaises(ValueError):
            primer_class_str_to_enum("")

    def test_strand_char_to_primer_class_str(self):
        self.assertEqual(strand_char_to_primer_class_str("+"), "LEFT")
        self.assertEqual(strand_char_to_primer_class_str("-"), "RIGHT")
        # Check unexpected
        with self.assertRaises(ValueError):
            strand_char_to_primer_class_str("")


class TestHeader(unittest.TestCase):
    def test_parse_headers_to_dict_simple(self):
        headers, _bedlines = read_bedfile(TEST_ATTRIBUTES_BEDFILE)

        attr_dict = parse_headers_to_dict(headers)
        self.assertDictEqual(
            attr_dict,
            {"MN908947.3": "sars-cov-2", "examplescheme": None, "gc": "fractiongc"},
        )

    def test_parse_headers_to_dict_probe(self):
        headers, _bedlines = read_bedfile(TEST_PROBE_BEDFILE)

        attr_dict = parse_headers_to_dict(headers)
        self.assertDictEqual(
            attr_dict,
            {
                "/3BHQ_1/": "BlackHoleQuencher1",
                "/56-FAM/": "FAM",
                "/5HEX/": "HEX",
                "example multiplexed-qPCR assay": None,
            },
        )


class TestAttributesFuncs(unittest.TestCase):
    def test_parse_primer_attributes_str_valid(self):
        # Valid str simple
        primer_attr = "pw=1"
        result = parse_primer_attributes_str(primer_attr)
        assert result is not None
        self.assertDictEqual(result, {PRIMER_WEIGHT_KEY: "1"})

        # valid multi
        primer_attr = "pw=1;ps=100;s=1"
        result = parse_primer_attributes_str(primer_attr)
        assert result is not None
        self.assertDictEqual(result, {"s": "1", PRIMER_WEIGHT_KEY: "1", "ps": "100"})

        # semi valid multi
        primer_attr = ";pw=1;ps=100;s=1;"
        result = parse_primer_attributes_str(primer_attr)
        assert result is not None
        self.assertDictEqual(result, {"s": "1", PRIMER_WEIGHT_KEY: "1", "ps": "100"})

        # Test whitespace + semi valid
        primer_attr = "; pw= 1;p s=100;s=1;"
        result = parse_primer_attributes_str(primer_attr)
        assert result is not None
        self.assertDictEqual(result, {"s": "1", PRIMER_WEIGHT_KEY: "1", "ps": "100"})

        # Test legacy weight
        primer_attr = "1.0"
        result = parse_primer_attributes_str(primer_attr)
        assert result is not None
        self.assertDictEqual(result, {PRIMER_WEIGHT_KEY: "1.0"})

    def test_parse_primer_attributes_str_invalid(self):
        # error ;
        primer_attr = ";"
        result = parse_primer_attributes_str(primer_attr)
        self.assertIsNone(result)

        # empty
        primer_attr = " "
        result = parse_primer_attributes_str(primer_attr)
        self.assertIsNone(result)

        # error malformed
        primer_attr = "pw="
        with self.assertRaises(ValueError) as context:
            parse_primer_attributes_str(primer_attr)
        self.assertIn("Malformed k=v pair:", str(context.exception))

        # error wrong form
        primer_attr = "pw1"
        with self.assertRaises(ValueError) as context:
            parse_primer_attributes_str(primer_attr)
        self.assertIn("Invalid PrimerAttributes:", str(context.exception))

        # error duplicate keys
        primer_attr = "pw=1;pw=2"
        with self.assertRaises(ValueError) as context:
            parse_primer_attributes_str(primer_attr)
        self.assertIn("Duplicate PrimerAttributes Key:", str(context.exception))

    def test_create_primer_attributes_str(self):
        # None case
        self.assertIsNone(create_primer_attributes_str(None))
        # empty dict
        self.assertIsNone(create_primer_attributes_str({}))
        # test valid single
        attr_str = create_primer_attributes_str({PRIMER_WEIGHT_KEY: "1"})
        self.assertEqual(attr_str, "pw=1")
        # test valid multi
        attr_str = create_primer_attributes_str({PRIMER_WEIGHT_KEY: "1", "ps": "As"})
        self.assertEqual(attr_str, "pw=1;ps=As")
        # test valid + sort order
        attr_str = create_primer_attributes_str(
            {"a": 100, PRIMER_WEIGHT_KEY: "1", "ps": "As"}
        )
        self.assertEqual(attr_str, "a=100;pw=1;ps=As")

        # Test strip whitespace
        attr_str = create_primer_attributes_str(
            {"a": 100, PRIMER_WEIGHT_KEY: "1", "ps": "A   s"}
        )
        self.assertEqual(attr_str, "a=100;pw=1;ps=As")

    def test_attributes_roundtrip(self):
        # Simple
        primer_attr = "pw=1"
        self.assertEqual(
            create_primer_attributes_str(parse_primer_attributes_str(primer_attr)),
            primer_attr,
        )
        # multi
        primer_attr = "a=100;pw=1;ps=As"
        self.assertEqual(
            create_primer_attributes_str(parse_primer_attributes_str(primer_attr)),
            primer_attr,
        )


class TestBedLine(unittest.TestCase):
    def setUp(self) -> None:
        self.bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        self.bedlinev2 = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT_1",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        return super().setUp()

    def test_create_bedline_with_strand_diff(self):
        """
        Test that
        """
        with self.assertRaises(ValueError):
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="scheme_1_LEFT",
                pool=1,
                strand="-",
                sequence="ACGT",
            )

    def test_bedline_create_left(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        # Provides values
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertIsNone(bedline.weight)

        # Derived values
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 1)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(
            bedline.amplicon_name,
            f"{bedline.amplicon_prefix}_{bedline.amplicon_number}",
        )
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertIsNone(bedline.primer_suffix)
        self.assertEqual(Strand.FORWARD, bedline.strand_class)

        self.assertEqual(bedline.ipool, 0)
        self.assertEqual(bedline.primer_class, PrimerClass.LEFT)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n",
        )

    def test_bedline_create_right(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_RIGHT",
            pool=1,
            strand="-",
            sequence="ACGT",
        )
        # Provides values
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_RIGHT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "-")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertIsNone(bedline.weight)

        # Derived values
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 1)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(
            bedline.amplicon_name,
            f"{bedline.amplicon_prefix}_{bedline.amplicon_number}",
        )
        self.assertIsNone(bedline.primer_suffix)
        self.assertEqual(Strand.REVERSE, bedline.strand_class)

        self.assertEqual(bedline.ipool, 0)
        self.assertEqual(bedline.primer_class, PrimerClass.RIGHT)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_RIGHT\t1\t-\tACGT\n",
        )

    def test_bedline_create_empty_weight(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
            attributes="",
        )
        # Provides values
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertIsNone(bedline.weight)

        # Derived values
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 1)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(bedline.amplicon_name, "scheme_1")
        self.assertEqual(bedline.ipool, 0)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n",
        )

    def test_bedline_create_probe(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_114_PROBE_2",
            pool=10,
            strand="+",
            sequence="ACGT",
            attributes="pw=10",
        )
        # Provides values
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_114_PROBE_2")
        self.assertEqual(bedline.pool, 10)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertEqual(bedline.weight, 10)
        assert bedline.attributes is not None
        self.assertDictEqual(bedline.attributes, {PRIMER_WEIGHT_KEY: 10})

        # Derived values
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 114)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(bedline.amplicon_name, "scheme_114")
        self.assertEqual(bedline.ipool, 9)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_114_PROBE_2\t10\t+\tACGT\tpw=10.0\n",
        )

    def test_bedline_create_weight(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
            attributes={PRIMER_WEIGHT_KEY: 1},
        )
        # Provides values
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertEqual(bedline.weight, 1.0)

        # Derived values
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 1)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(bedline.ipool, 0)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\tpw=1.0\n",
        )

    def test_bedline_create_attr(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
            attributes="pw=1;ps=2;dfa=1000",
        )
        # Provides values
        self.assertEqual(bedline.chrom, "chr1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.primername, "scheme_1_LEFT")
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.sequence, "ACGT")
        self.assertEqual(bedline.weight, 1.0)
        assert bedline.attributes is not None

        # Check attr
        self.assertEqual(bedline.attributes["pw"], 1)
        self.assertEqual(bedline.attributes["ps"], "2")
        self.assertEqual(bedline.attributes["dfa"], "1000")

        # Derived values
        self.assertEqual(bedline.length, 100)
        self.assertEqual(bedline.amplicon_number, 1)
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(bedline.ipool, 0)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\tpw=1.0;ps=2;dfa=1000\n",
        )

        # Change an attr
        bedline.attributes["test"] = 100
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\tpw=1.0;ps=2;dfa=1000;test=100\n",
        )

        # Remove pw
        bedline.attributes.pop("pw")
        self.assertIsNone(bedline.weight)
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\tps=2;dfa=1000;test=100\n",
        )

    def test_create_invalid_bedline(self):
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
        # Invalid weight should raise ValueError
        with self.assertRaises(ValueError):
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="scheme_1_LEFT",
                pool=1,
                strand="+",
                sequence="ACGT",
                attributes={PRIMER_WEIGHT_KEY: -1.0},
            )
        # str weight should raise ValueError
        with self.assertRaises(ValueError):
            BedLine(
                chrom="chr1",
                start=100,
                end=200,
                primername="scheme_1_LEFT",
                pool=1,
                strand="+",
                sequence="ACGT",
                attributes={PRIMER_WEIGHT_KEY: "A"},
            )

    def test_bedline_parse_params(self):
        bedline = BedLine(
            chrom=1,  # int chrom # type: ignore
            start="100",  # str start
            end="200",  # str end
            primername="scheme_1_LEFT",
            pool="1",  # str pool
            strand="+",
            sequence="atcg",  # lowercase sequence
        )
        self.assertEqual(bedline.chrom, "1")
        self.assertEqual(bedline.start, 100)
        self.assertEqual(bedline.end, 200)
        self.assertEqual(bedline.pool, 1)
        self.assertEqual(bedline.sequence, "ATCG")

    def test_bedline_parse_params_invalid(self):
        valid_bedline = BedLine(
            chrom="chr1",
            start="100",
            end="200",
            primername="scheme_1_LEFT",
            pool="1",
            strand="+",
            sequence="ATCG",
        )

        # Invalid pool should raise ValueError
        with self.assertRaises(ValueError):
            valid_bedline.pool = "0"

        # Invalid strand should raise ValueError
        with self.assertRaises(ValueError):
            valid_bedline.strand = "invalid"

        # Invalid start should raise ValueError
        with self.assertRaises(ValueError):
            valid_bedline.start = "invalid"

        # Invalid end should raise ValueError
        with self.assertRaises(ValueError):
            valid_bedline.end = "invalid"

        # Invalid primername should raise ValueError
        with self.assertRaises(ValueError):
            valid_bedline.primername = "invalid"

    def test_bedline_change_class(self):
        bl = self.bedlinev2

        # CHeck is left
        self.assertEqual(bl.primername, "scheme_1_LEFT_1")
        self.assertEqual(bl.strand, Strand.FORWARD.value)

        # Check LEFT -> PROBE
        bl.primer_class = PrimerClass.PROBE.value
        self.assertEqual(bl.primername, "scheme_1_PROBE_1")
        self.assertEqual(bl.strand, Strand.FORWARD.value)  # Check strand unchanged

        # Check PROBE -> RIGHT
        with self.assertRaises(ValueError):
            bl.primer_class = PrimerClass.RIGHT.value
        # Change the strand
        bl.strand = Strand.REVERSE.value

        # PROBE -> RIGHT
        bl.primer_class = PrimerClass.RIGHT.value
        self.assertEqual(bl.primername, "scheme_1_RIGHT_1")
        self.assertEqual(bl.strand, Strand.REVERSE.value)  # Check strand unchanged

        # RIGHT -> LEFT
        with self.assertRaises(ValueError):
            bl.primer_class = PrimerClass.LEFT.value
        # Try force change
        bl.force_change(PrimerClass.LEFT.value, Strand.FORWARD.value)
        self.assertEqual(bl.primername, "scheme_1_LEFT_1")
        self.assertEqual(bl.strand, Strand.FORWARD.value)

        # Force change LEFT -> RIGHT:
        with self.assertRaises(ValueError):
            bl.primer_class = PrimerClass.RIGHT.value
        # Try force change
        bl.force_change(PrimerClass.RIGHT.value, Strand.REVERSE.value)
        self.assertEqual(bl.primername, "scheme_1_RIGHT_1")
        self.assertEqual(bl.strand, Strand.REVERSE.value)

    def test_to_bed(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n",
        )
        # Provide weight
        bedline.weight = 1.0
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\tpw=1.0\n",
        )

    def test_to_bed_probe(self):
        bedline = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_PROBE",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_PROBE\t1\t+\tACGT\n",
        )
        # Provide weight
        bedline.weight = 1.0
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_PROBE\t1\t+\tACGT\tpw=1.0\n",
        )

        bedline = update_primernames([bedline])[0]
        self.assertEqual(
            bedline.to_bed(),
            "chr1\t100\t200\tscheme_1_PROBE_1\t1\t+\tACGT\tpw=1.0\n",
        )

    def test_chrom_set(self):
        bedline = self.bedline
        # Object
        with self.assertRaises(ValueError):
            bedline.chrom = []

        # Valid
        bedline.chrom = "test"
        self.assertEqual(bedline.chrom, "test")

        # Strip
        bedline.chrom = "test1    "
        self.assertEqual(bedline.chrom, "test1")

        # Internal whitespace
        bedline.chrom = "test  2"
        self.assertEqual(bedline.chrom, "test2")

        # Int convert
        bedline.chrom = 1
        self.assertEqual(bedline.chrom, "1")

    def test_start_setter(self):
        # Test setting a valid positive integer
        bedline = self.bedline
        bedline.start = 50
        self.assertEqual(bedline.start, 50)

        # Test setting zero
        bedline.start = 0
        self.assertEqual(bedline.start, 0)

        # Str convert
        bedline.start = "75"
        self.assertEqual(bedline.start, 75)

        # Float convert
        bedline.start = 80.5
        self.assertEqual(bedline.start, 80)

        # Negative
        with self.assertRaises(ValueError) as context:
            bedline.start = -5
        self.assertIn(
            "start must be greater than or equal to 0", str(context.exception)
        )

        # Negative str
        with self.assertRaises(ValueError) as context:
            bedline.start = "-5"
        self.assertIn(
            "start must be greater than or equal to 0", str(context.exception)
        )

        # Object
        with self.assertRaises(ValueError) as context:
            bedline.start = [1, 2, 3]
        self.assertIn("start must be an int", str(context.exception))

    def test_end_setter(self):
        # Test setting a valid positive integer
        bedline = self.bedline
        bedline.end = 250
        self.assertEqual(bedline.end, 250)

        # Test setting to start value
        bedline.end = 100
        self.assertEqual(bedline.end, 100)

        # Str convert
        bedline.end = "175"
        self.assertEqual(bedline.end, 175)

        # Float convert
        bedline.end = 180.5
        self.assertEqual(bedline.end, 180)

        # Negative
        with self.assertRaises(ValueError) as context:
            bedline.end = -5
        self.assertIn("end must be greater than or equal to 0", str(context.exception))

        # Negative str
        with self.assertRaises(ValueError) as context:
            bedline.end = "-5"
        self.assertIn("end must be greater than or equal to 0", str(context.exception))

        # Object
        with self.assertRaises(ValueError) as context:
            bedline.end = [1, 2, 3]
        self.assertIn("end must be an int", str(context.exception))

    def test_weight_setter(self):
        bedline = self.bedline
        # Str non int
        with self.assertRaises(ValueError):
            bedline.weight = "!"
        # Negative
        with self.assertRaises(ValueError):
            bedline.weight = -1.0
        # Negative Str
        with self.assertRaises(ValueError):
            bedline.weight = "-1.0"

        # Valid
        bedline.weight = 1.0
        self.assertEqual(bedline.weight, 1.0)

        # Explicit None
        bedline.weight = None
        self.assertIsNone(bedline.weight)

        # Parse str
        bedline.weight = "2.0"
        self.assertEqual(bedline.weight, 2.0)

        # empty str -> None
        bedline.weight = ""
        self.assertIsNone(bedline.weight)

    def test_update_primername_strand(self):
        bedline = self.bedline
        bedline.primername = "scheme_10_RIGHT_2"

        # Check all updated
        self.assertEqual(bedline.amplicon_prefix, "scheme")
        self.assertEqual(bedline.amplicon_number, 10)
        self.assertEqual(bedline.strand, "-")
        self.assertEqual(bedline.primer_class, PrimerClass.RIGHT)

    def test_pool_setter(self):
        bedline = self.bedline

        # str
        with self.assertRaises(ValueError) as cm:
            bedline.pool = "A"
        self.assertIn("pool must be an int", str(cm.exception))

        # list
        with self.assertRaises(ValueError) as cm:
            bedline.pool = []
        self.assertIn("pool must be an int", str(cm.exception))

    def test_primername_suffix_setter(self):
        pass

    def test_primername_setter(self):
        # Create a basic BedLine instance
        bedline = self.bedline

        # Test V1 format primername
        bedline.primername = "test_2_LEFT"
        self.assertEqual(bedline.amplicon_prefix, "test")
        self.assertEqual(bedline.amplicon_number, 2)
        self.assertEqual(bedline.strand, "+")
        self.assertIsNone(bedline.primer_suffix)
        self.assertEqual(bedline.primername, "test_2_LEFT")

        # Test V1 format with alt suffix
        bedline.primername = "test_3_LEFT_alt1"
        self.assertEqual(bedline.amplicon_prefix, "test")
        self.assertEqual(bedline.amplicon_number, 3)
        self.assertEqual(bedline.strand, "+")
        self.assertEqual(bedline.primer_suffix, "alt1")
        self.assertEqual(bedline.primername, "test_3_LEFT_alt1")

        # Test V2 format
        bedline.primername = "test_4_RIGHT_2"
        self.assertEqual(bedline.amplicon_prefix, "test")
        self.assertEqual(bedline.amplicon_number, 4)
        self.assertEqual(bedline.strand, "-")
        self.assertEqual(bedline.primer_suffix, 2)
        self.assertEqual(bedline.primername, "test_4_RIGHT_2")

        # Test with hyphenated amplicon prefix
        bedline.primername = "SARS-CoV-2_5_LEFT"
        self.assertEqual(bedline.amplicon_prefix, "SARS-CoV-2")
        self.assertEqual(bedline.amplicon_number, 5)
        self.assertEqual(bedline.strand, "+")
        self.assertIsNone(bedline.primer_suffix)
        self.assertEqual(bedline.primername, "SARS-CoV-2_5_LEFT")

        # Invalid primername (wrong format)
        with self.assertRaises(ValueError) as context:
            bedline.primername = "invalid_primername"
        self.assertIn("Invalid primername", str(context.exception))

        # Invalid primername (missing LEFT/RIGHT)
        with self.assertRaises(ValueError) as context:
            bedline.primername = "test_6_MIDDLE"
        self.assertIn("Invalid primername", str(context.exception))

        # Invalid primername (non-numeric amplicon number)
        with self.assertRaises(ValueError) as context:
            bedline.primername = "test_A_LEFT"
        self.assertIn("Invalid primername", str(context.exception))

    def test_update_primername(self):
        bedline = self.bedline
        # Update amplicon number
        bedline.amplicon_number = 10
        self.assertEqual(bedline.primername, "scheme_10_LEFT")
        # Update amplicon prefix
        bedline.amplicon_prefix = "test"
        self.assertEqual(bedline.primername, "test_10_LEFT")
        # Update the primer suffix
        bedline.primer_suffix = "alt1"
        self.assertEqual(bedline.primername, "test_10_LEFT_alt1")

        # Test invalid primer suffix
        with self.assertRaises(ValueError):
            bedline.primer_suffix = "test"
        with self.assertRaises(ValueError):
            bedline.primer_suffix = " "
        with self.assertRaises(ValueError):
            bedline.primer_suffix = "none"

        # Test prefix
        with self.assertRaises(ValueError):
            bedline.amplicon_prefix = "10_left"
        with self.assertRaises(ValueError):
            bedline.amplicon_prefix = "."

        # test amplicon number
        with self.assertRaises(ValueError):
            bedline.amplicon_number = "A"

    def test_primername_version(self):
        bedline = self.bedline
        # Check old
        self.assertEqual(bedline.primername_version, PrimerNameVersion.V1)

        # Update
        bedline.primername = "scheme_1_LEFT_1"
        self.assertEqual(bedline.primername_version, PrimerNameVersion.V2)

    def test_primer_suffix_setter(self):
        bedline = self.bedline

        # Invalid primer_suffix (wrong object)
        with self.assertRaises(ValueError) as context:
            bedline.primer_suffix = []  # type: ignore
        self.assertIn("Invalid primer_suffix", str(context.exception))

        # Invalid primer_suffix (wrong object)
        with self.assertRaises(ValueError) as context:
            bedline.primer_suffix = -1
        self.assertIn(
            "primer_suffix must be greater than or equal to 0", str(context.exception)
        )

        # Invalid v1 primer_suffix
        with self.assertRaises(ValueError) as context:
            bedline.primer_suffix = "A"
        self.assertIn("Invalid V1 primer_suffix:", str(context.exception))

        # Set None
        bedline.primer_suffix = None
        self.assertIsNone(bedline.primer_suffix)

    def test_attribute_setter(self):
        bedline = self.bedline

        # Set nonsense
        with self.assertRaises(ValueError) as context:
            bedline.attributes = []  # type: ignore
        self.assertIn("Invalid primer attributes.", str(context.exception))

        # Set empty dict
        bedline.attributes = {}
        self.assertEqual(bedline.attributes, {})
        # Ensure empty dict is not written to bed
        self.assertEqual(
            bedline.to_bed(), "chr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )

        # Set string. Test pw is converted to float
        bedline.attributes = "pw=1;ps=2"
        assert bedline.attributes is not None
        self.assertDictEqual(bedline.attributes, {"pw": 1.0, "ps": "2"})  # type: ignore

        # Set an invalid pw
        with self.assertRaises(ValueError) as context:
            bedline.attributes = "pw=A;"
        self.assertIn("weight must be a float", str(context.exception))

        # Check attribute will be empty dict
        bedline.attributes = None
        self.assertEqual(bedline.attributes, {})

        # Check the (if bedline.attributes:) pattern works
        self.assertFalse(bedline.attributes)


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

    def test_create_bedline_invalid(self):
        with self.assertRaises(IndexError):
            create_bedline([""])


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

    def test_read_weight_bedline(self):
        headers, bedlines = read_bedfile(TEST_WEIGHTS_BEDFILE)
        # Check for empty headers
        self.assertEqual(headers, [])

        # Check bedlines have weights
        for bedline in bedlines:
            self.assertIsNotNone(bedline.weight)

    def test_read_attr_bedline(self):
        headers, bedlines = read_bedfile(TEST_ATTRIBUTES_BEDFILE)
        # Check for 3 headers
        self.assertEqual(len(headers), 3)

        # Check bedlines have weights and attrs
        for bedline in bedlines:
            self.assertIsNotNone(bedline.weight)
            self.assertIsNotNone(bedline.attributes)
            assert bedline.attributes is not None
            # Check for two attrs
            self.assertEqual(len(bedline.attributes), 2)


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
        # Write non-weighted bedline
        write_bedfile(self.output_bed_path, ["#header1"], [bedline])
        with open(self.output_bed_path) as f:
            content = f.read()
        self.assertEqual(
            content, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\n"
        )
        # Write weighted bedline
        bedline.weight = 1.0
        write_bedfile(self.output_bed_path, ["#header1"], [bedline])
        with open(self.output_bed_path) as f:
            content = f.read()
        self.assertEqual(
            content, "#header1\nchr1\t100\t200\tscheme_1_LEFT\t1\t+\tACGT\tpw=1.0\n"
        )

    def tearDown(self) -> None:
        # Remove output file
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


class TestGroupByPool(unittest.TestCase):
    def test_group_by_pool(self):
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

        grouped = group_by_pool([bedline1, bedline2])
        self.assertEqual(len(grouped), 2)

        # Check for correct pools
        self.assertEqual(grouped.keys(), {1, 2})

        # Change bedline1 pool to 2
        bedline1.pool = 20

        grouped = group_by_pool([bedline1, bedline2])
        self.assertEqual(grouped.keys(), {2, 20})


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


class TestGroupByDirection(unittest.TestCase):
    def test_group_by_direction(self):
        bedline1 = BedLine(
            chrom="chr1",
            start=100,
            end=200,
            primername="scheme_1_LEFT_1",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        bedline2 = BedLine(
            chrom="chr1",
            start=150,
            end=250,
            primername="scheme_1_RIGHT_1",
            pool=2,
            strand="-",
            sequence="ACGT",
        )
        bedline3 = BedLine(
            chrom="chr1",
            start=175,
            end=225,
            primername="scheme_1_PROBE_1",
            pool=1,
            strand="+",
            sequence="ACGT",
        )
        grouped = group_by_class([bedline1, bedline2, bedline3])
        self.assertEqual(len(grouped), 3)
        self.assertEqual(len(grouped[PrimerClass.LEFT.value]), 1)
        self.assertEqual(len(grouped[PrimerClass.RIGHT.value]), 1)
        self.assertEqual(len(grouped[PrimerClass.PROBE.value]), 1)
        self.assertEqual(grouped[PrimerClass.LEFT.value], [bedline1])
        self.assertEqual(grouped[PrimerClass.RIGHT.value], [bedline2])
        self.assertEqual(grouped[PrimerClass.PROBE.value], [bedline3])

    def test_group_by_direction_file(self):
        headers, bedlines = read_bedfile(TEST_BEDFILE)
        grouped = group_by_class(bedlines)
        self.assertEqual(len(grouped), 2)
        self.assertEqual(len(grouped[PrimerClass.LEFT.value]), 3)
        self.assertEqual(len(grouped[PrimerClass.RIGHT.value]), 3)

        # Check that all LEFT primers have correct direction
        for bedline in grouped[PrimerClass.LEFT.value]:
            self.assertEqual(bedline.direction_str, "LEFT")

        # Check that all RIGHT primers have correct direction
        for bedline in grouped[PrimerClass.RIGHT.value]:
            self.assertEqual(bedline.direction_str, "RIGHT")


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

    def test_merge_primers_single(self):
        bedlines = [
            BedLine(
                chrom="chr1",
                start=100,
                end=120,
                primername="test_1_RIGHT_1",
                pool=1,
                strand="-",
                sequence="ACGT",
            ),
            BedLine(
                chrom="chr1",
                start=110,
                end=130,
                primername="test_1_RIGHT_2",
                pool=1,
                strand="-",
                sequence="ACGT",
            ),
        ]
        merged_bedlines = merge_primers(bedlines)
        # Check merged bedline
        self.assertEqual(len(merged_bedlines), 1)

        # Check merged bedline attributes
        merged_bedline = merged_bedlines[0]
        self.assertEqual(merged_bedline.chrom, "chr1")
        self.assertEqual(merged_bedline.start, 100)
        self.assertEqual(merged_bedline.end, 130)
        self.assertEqual(merged_bedline.primername, "test_1_RIGHT_1")
        self.assertEqual(merged_bedline.pool, 1)
        self.assertEqual(merged_bedline.strand, "-")

    def test_merge_primers_nothing(self):
        bedlines = [
            BedLine(
                chrom="chr1",
                start=100,
                end=120,
                primername="test_1_LEFT_1",
                pool=1,
                strand="+",
                sequence="ACGT",
            ),
            BedLine(
                chrom="chr1",
                start=110,
                end=130,
                primername="test_1_RIGHT_2",
                pool=1,
                strand="-",
                sequence="ACGT",
            ),
        ]
        merged_bedlines = merge_primers(bedlines)

        self.assertEqual(len(merged_bedlines), 2)

    def test_merge_primers_empty(self):
        merged_bedlines = merge_primers([])
        self.assertEqual(len(merged_bedlines), 0)

    def test_expand_bedlines(self):
        bedlines = [
            BedLine(
                chrom="chr1",
                start=100,
                end=120,
                primername="test_1_LEFT_1",
                pool=1,
                strand="+",
                sequence="ACNT",
            ),
            BedLine(
                chrom="chr1",
                start=110,
                end=130,
                primername="test_1_RIGHT_1",
                pool=1,
                strand="-",
                sequence="AYGT",
            ),
            BedLine(
                chrom="chr1",
                start=110,
                end=130,
                primername="test_2_LEFT_1",
                pool=1,
                strand="+",
                sequence="ACTGATCGAC",
            ),
        ]
        expected_names = [
            "test_1_LEFT_1",
            "test_1_LEFT_2",
            "test_1_LEFT_3",
            "test_1_LEFT_4",
            "test_1_RIGHT_1",
            "test_1_RIGHT_2",
            "test_2_LEFT_1",
        ]
        bedlines = expand_bedlines(bedlines)

        primer_names = [bl.primername for bl in bedlines]

        self.assertEqual(expected_names, primer_names)


if __name__ == "__main__":
    unittest.main()
