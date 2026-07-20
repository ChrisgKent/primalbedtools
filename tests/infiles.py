import pathlib

# Bedfiles
TEST_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.bed"
TEST_V2_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.v2.bed"
TEST_WEIGHTS_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.weights.bed"
TEST_WEIGHTS_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.weights.bed"
TEST_ATTRIBUTES_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.attributes.bed"
TEST_PROBE_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.probe.bed"
TEST_PANEL_BEDFILE = pathlib.Path(__file__).parent / "inputs/panel.input.bed"
TEST_MIXED_PREFIX_BEDFILE = (
    pathlib.Path(__file__).parent / "inputs/test.mixed_prefix.bed"
)

# fasta
FASTA_PATH = pathlib.Path(__file__).parent / "inputs/msa.input.fasta"
