import contextlib
import io
import sys
import unittest
from unittest import mock

from primalbedtools.logs import reset_cli_logging
from primalbedtools.main import main
from tests.infiles import TEST_BEDFILE, TEST_MIXED_PREFIX_BEDFILE


class CliTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        # Stop the stderr handler leaking into other tests in this process
        reset_cli_logging()
        return super().tearDown()

    def run_cli(self, argv):
        """Run main() with argv, returning (exit_code, stdout, stderr)."""
        out, err = io.StringIO(), io.StringIO()
        with mock.patch.object(sys, "argv", ["pbt"] + argv):
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                with self.assertRaises(SystemExit) as cm:
                    main()
        return cm.exception.code, out.getvalue(), err.getvalue()


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


if __name__ == "__main__":
    unittest.main()
