import contextlib
import io
import os
import sys
import tempfile
import unittest
from unittest import mock

from primalbedtools.logs import reset_cli_logging
from primalbedtools.main import main
from tests.infiles import (
    REFERENCE_PATH,
    TEST_ATTRIBUTES_BEDFILE,
    TEST_BEDFILE,
    TEST_MIXED_PREFIX_BEDFILE,
    TEST_PRIMER_BEDFILE,
)


class CliTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        # Stop the stderr handler leaking into other tests in this process
        reset_cli_logging()
        return super().tearDown()

    def run_cli(self, argv, stdin=None):
        """Run main() with argv, returning (exit_code, stdout, stderr)."""
        out, err = io.StringIO(), io.StringIO()
        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(sys, "argv", ["pbt"] + argv))
            stack.enter_context(contextlib.redirect_stdout(out))
            stack.enter_context(contextlib.redirect_stderr(err))
            if stdin is not None:
                stack.enter_context(mock.patch.object(sys, "stdin", io.StringIO(stdin)))
            # Most subcommands exit(0); returning normally is also a success
            code = 0
            try:
                main()
            except SystemExit as e:
                code = e.code
        return code, out.getvalue(), err.getvalue()


class TestAmpliconCli(CliTestCase):
    def test_stdout_is_valid_bed(self):
        # The prefix warning must not be interleaved with the BED records
        code, out, err = self.run_cli(["amplicon", str(TEST_MIXED_PREFIX_BEDFILE)])

        self.assertEqual(code, 0)
        self.assertTrue(out)
        for line in out.splitlines():
            self.assertEqual(len(line.split("\t")), 5, f"not a BED record: {line!r}")
        self.assertIn("more than one prefix", err)

    def test_joined_prefix_in_amplicon_name(self):
        _code, out, _err = self.run_cli(["amplicon", str(TEST_MIXED_PREFIX_BEDFILE)])

        names = [line.split("\t")[3] for line in out.splitlines()]
        self.assertEqual(
            names,
            [
                "SARS-CoV-2-F-SARS-CoV-2-R_1",
                "SARS-CoV-2-F-SARS-CoV-2-R_2",
                "SARS-CoV-2-F-SARS-CoV-2-R_3",
            ],
        )

    def test_quiet_suppresses_warning(self):
        _code, out, err = self.run_cli(
            ["-q", "amplicon", str(TEST_MIXED_PREFIX_BEDFILE)]
        )

        self.assertEqual(err, "")
        self.assertTrue(out)

    def test_verbose_emits_per_amplicon_debug(self):
        _code, _out, err = self.run_cli(
            ["-v", "amplicon", str(TEST_MIXED_PREFIX_BEDFILE)]
        )

        self.assertIn("DEBUG", err)
        self.assertIn("differing prefixes", err)

    def test_clean_bedfile_is_silent(self):
        _code, out, err = self.run_cli(["amplicon", str(TEST_BEDFILE)])

        self.assertEqual(err, "")
        self.assertIn("SARS-CoV-2_1", out)


class TestValidateCli(CliTestCase):
    def test_clean_bedfile_emits_no_prefix_noise(self):
        _code, _out, err = self.run_cli(["validate_bedfile", str(TEST_BEDFILE)])

        self.assertNotIn("prefix", err)


class TestFromCsvCli(CliTestCase):
    def test_round_trips_via_csv(self):
        _code, csv_out, _err = self.run_cli(["csv", str(TEST_ATTRIBUTES_BEDFILE)])

        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as f:
            f.write(csv_out)
            csv_path = f.name
        self.addCleanup(os.unlink, csv_path)

        code, bed_out, _err = self.run_cli(["from-csv", csv_path])

        self.assertEqual(code, 0)
        # Headers are dropped, so compare against the bedlines only
        expected = "".join(
            line + "\n"
            for line in TEST_ATTRIBUTES_BEDFILE.read_text().splitlines()
            if not line.startswith("#")
        )
        self.assertEqual(bed_out, expected)


class TestSubcommandSmoke(CliTestCase):
    """Every subcommand runs, exits 0, and writes only its data to stdout."""

    def assert_clean_run(self, argv, expected_in_stdout):
        code, out, err = self.run_cli(argv)

        self.assertEqual(code, 0, f"{argv} exited {code}: {err}")
        self.assertIn(expected_in_stdout, out)
        self.assertEqual(err, "", f"{argv} wrote to stderr: {err}")

    def test_sort(self):
        self.assert_clean_run(["sort", str(TEST_BEDFILE)], "SARS-CoV-2_1_LEFT_1")

    def test_update(self):
        self.assert_clean_run(["update", str(TEST_BEDFILE)], "SARS-CoV-2_1_LEFT_1")

    def test_merge(self):
        self.assert_clean_run(["merge", str(TEST_BEDFILE)], "SARS-CoV-2_1_LEFT")

    def test_format(self):
        self.assert_clean_run(["format", str(TEST_BEDFILE)], "SARS-CoV-2_1_LEFT_1")

    def test_fasta(self):
        self.assert_clean_run(["fasta", str(TEST_BEDFILE)], ">SARS-CoV-2_1_LEFT_1")

    def test_csv(self):
        self.assert_clean_run(["csv", str(TEST_BEDFILE)], "chrom,start,end")

    def test_downgrade(self):
        self.assert_clean_run(["downgrade", str(TEST_BEDFILE)], "SARS-CoV-2_1_LEFT")

    def test_amplicon_primertrim(self):
        self.assert_clean_run(["amplicon", "-t", str(TEST_BEDFILE)], "SARS-CoV-2_1")

    def test_validate_bedfile(self):
        code, _out, err = self.run_cli(["validate_bedfile", str(TEST_BEDFILE)])

        self.assertEqual(code, 0)
        self.assertEqual(err, "")


class TestDiffCli(CliTestCase):
    def test_both_inputs_from_stdin_is_rejected(self):
        # stdin can only be read once, so this silently diffed against nothing
        code, _out, err = self.run_cli(
            ["diff", "-", "-"], stdin=TEST_BEDFILE.read_text()
        )

        self.assertNotEqual(code, 0)
        self.assertIn("only one of bedfile1/bedfile2", err)

    def test_one_input_from_stdin_works(self):
        code, out, _err = self.run_cli(
            ["diff", "-", str(TEST_BEDFILE)], stdin=TEST_BEDFILE.read_text()
        )

        self.assertEqual(code, 0)
        self.assertEqual(out, "")


class TestStdinInput(CliTestCase):
    """A bed path of "-" reads from stdin."""

    def setUp(self) -> None:
        self.bed_text = TEST_BEDFILE.read_text()
        return super().setUp()

    def test_format_round_trips_through_stdin(self):
        _code, from_stdin, _err = self.run_cli(["format", "-"], stdin=self.bed_text)
        _code, from_file, _err = self.run_cli(["format", str(TEST_BEDFILE)])

        self.assertEqual(from_stdin, from_file)
        self.assertIn("SARS-CoV-2_1_LEFT_1", from_stdin)

    def test_amplicon_reads_stdin(self):
        _code, out, _err = self.run_cli(["amplicon", "-"], stdin=self.bed_text)

        self.assertIn("SARS-CoV-2_1", out)
        for line in out.splitlines():
            self.assertEqual(len(line.split("\t")), 5)

    def test_sort_reads_stdin(self):
        _code, out, _err = self.run_cli(["sort", "-"], stdin=self.bed_text)

        self.assertIn("SARS-CoV-2_1_LEFT_1", out)

    def test_validate_reads_stdin_once(self):
        # validate used to re-open the path after main() had already read it,
        # which consumes stdin twice and finds nothing the second time.
        code, _out, _err = self.run_cli(
            ["validate", "-", str(REFERENCE_PATH)],
            stdin=TEST_PRIMER_BEDFILE.read_text(),
        )

        self.assertEqual(code, 0)

    def test_validate_from_file_still_works(self):
        code, _out, _err = self.run_cli(
            ["validate", str(TEST_PRIMER_BEDFILE), str(REFERENCE_PATH)]
        )

        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
