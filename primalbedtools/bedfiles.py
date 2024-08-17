import enum
import pathlib
import re

from pydantic import BaseModel, ConfigDict, field_validator

# Regular expressions for primer names
V2_PRIMERNAME = r"^[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)_[0-9]+$"
V1_PRIMERNAME = r"^[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)(_ALT[0-9]*|_alt[0-9]*)*$"


def check_primername(primername: str) -> bool:
    """
    Check if a primername is valid.
    """
    return bool(
        re.match(V1_PRIMERNAME, primername) or re.match(V2_PRIMERNAME, primername)
    )


class StrandEnum(enum.Enum):
    FORWARD = "+"
    REVERSE = "-"


class BedLine(BaseModel):
    """
    A BedLine object represents a single line in a BED file.

    Attributes:
    - chrom: str
    - start: int
    - end: int
    - primername: str
    - pool: int # 1-based pool number use ipool for 0-based pool number
    - strand: StrandEnum
    - sequence : str
    """

    # pydantic
    model_config = ConfigDict(use_enum_values=True, str_strip_whitespace=True)

    # properties
    chrom: str
    start: int
    end: int
    primername: str
    pool: int
    strand: StrandEnum | str
    sequence: str

    @field_validator("primername")
    @classmethod
    def check_primername(cls, v):
        if not check_primername(v):
            raise ValueError(f"Invalid primername: {v}")
        return v

    # calculated properties
    @property
    def length(self):
        return self.end - self.start

    @property
    def amplicon_number(self) -> int:
        return int(self.primername.split("_")[1])

    @property
    def amplicon_prefix(self) -> str:
        return self.primername.split("_")[0]

    @property
    def ipool(self) -> int:
        """Return the 0-based pool number"""
        return self.pool - 1

    def to_bed(self) -> str:
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.primername}\t{self.pool}\t{self.strand}\t{self.sequence}\n"


class BedLineParser:
    """
    Collection of methods for BED file IO.
    """

    @staticmethod
    def from_file(bedfile: str | pathlib.Path) -> tuple[list[str], list[BedLine]]:
        """
        Read and parse a BED file and return a tuple of headers and BedLine objects.
        : param bedfile: str | pathlib.Path
        : return: tuple[list[str], list[BedLine]]
        """
        return read_bedfile(bedfile=bedfile)

    @staticmethod
    def from_str(bedfile_str: str) -> tuple[list[str], list[BedLine]]:
        """
        Parse a BED string and return a tuple of headers and BedLine objects.
        : param bedfile_str: str
        : return: tuple[list[str], list[BedLine]]
        """
        bedfile_lines = bedfile_str.strip().split("\n")
        headers = []
        bedlines = []
        for line in bedfile_lines:
            line = line.strip()
            if line.startswith("#"):
                headers.append(line)
            elif line:
                bedlines.append(create_bedline(line.split("\t")))
        return headers, bedlines

    @staticmethod
    def to_str(headers: list[str] | None, bedlines: list[BedLine]) -> str:
        """
        Creates a BED string from the headers and BedLine objects.
        : param headers: list[str] | None
        : param bedlines: list[BedLine]
        : return: str
        """
        return create_bedfile_str(headers, bedlines)

    @staticmethod
    def to_file(
        bedfile: str | pathlib.Path, headers: list[str] | None, bedlines: list[BedLine]
    ) -> None:
        """
        Creates a BED file from the headers and BedLine objects.
        : param bedfile: str | pathlib.Path
        : param headers: list[str] | None
        : param bedlines: list[BedLine]
        """
        write_bedfile(bedfile, headers, bedlines)


def create_bedline(bedline: list[str]) -> BedLine:
    return BedLine(
        chrom=bedline[0],
        start=int(bedline[1]),
        end=int(bedline[2]),
        primername=bedline[3],
        pool=int(bedline[4]),
        strand=StrandEnum(bedline[5]),
        sequence=bedline[6],
    )


def read_bedfile(bedfile: str | pathlib.Path) -> tuple[list[str], list[BedLine]]:
    headers = []
    bedlines = []
    with open(bedfile) as f:
        for line in f.readlines():
            line = line.strip()

            if line.startswith("#"):
                headers.append(line)
                continue
            else:
                bedlines.append(create_bedline(line.split("\t")))

    return headers, bedlines


def create_bedfile_str(headers: list[str] | None, bedlines: list[BedLine]) -> str:
    bedfile_str: list[str] = []
    if headers:
        for header in headers:
            # add # if not present
            if not header.startswith("#"):
                header = "#" + header
            bedfile_str.append(header + "\n")
    # Add bedlines
    for bedline in bedlines:
        bedfile_str.append(bedline.to_bed())

    return "".join(bedfile_str)


def write_bedfile(
    bedfile: str | pathlib.Path, headers: list[str] | None, bedlines: list[BedLine]
):
    with open(bedfile, "w") as f:
        f.write(create_bedfile_str(headers, bedlines))


def group_by_chrom(list_bedlines: list[BedLine]) -> dict[str, list[BedLine]]:
    """
    Group a list of BedLine objects by chrom attribute.
    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.chrom not in bedlines_dict:
            bedlines_dict[bedline.chrom] = []
        bedlines_dict[bedline.chrom].append(bedline)
    return bedlines_dict


def group_by_amplicon_number(list_bedlines: list[BedLine]) -> dict[int, list[BedLine]]:
    """
    Group a list of BedLine objects by amplicon number.
    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.amplicon_number not in bedlines_dict:
            bedlines_dict[bedline.amplicon_number] = []
        bedlines_dict[bedline.amplicon_number].append(bedline)
    return bedlines_dict


def group_by_strand(
    list_bedlines: list[BedLine],
) -> dict[StrandEnum | str, list[BedLine]]:
    """
    Group a list of BedLine objects by strand.
    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.strand not in bedlines_dict:
            bedlines_dict[bedline.strand] = []
        bedlines_dict[bedline.strand].append(bedline)
    return bedlines_dict


def group_primer_pairs(bedlines: list[BedLine]) -> list[tuple[BedLine, BedLine]]:
    """
    Generate primer pairs from a list of BedLine objects.
    Groups by chrom, then by amplicon number, then pairs forward and reverse primers.
    """
    primer_pairs = []

    # Group by chrom
    for chrom_bedlines in group_by_chrom(bedlines).values():
        # Group by amplicon number
        for amplicon_number_bedlines in group_by_amplicon_number(
            chrom_bedlines
        ).values():
            # Generate primer pairs
            strand_to_bedlines = group_by_strand(amplicon_number_bedlines)
            primer_pairs.append(
                (
                    strand_to_bedlines.get(StrandEnum.FORWARD, []),
                    strand_to_bedlines.get(StrandEnum.REVERSE, []),
                )
            )

    return primer_pairs


class PrimerPair:
    """
    A PrimerPair object represents an amplicon with forward and reverse primers.
    """

    fbedlines: list[BedLine]
    rbedlines: list[BedLine]

    chrom: str
    pool: int
    amplicon_number: int
    prefix: str

    def __init__(self, fbedlines: list[BedLine], rbedlines: list[BedLine]):
        self.fbedlines = fbedlines
        self.rbedlines = rbedlines

        all_lines = fbedlines + rbedlines

        # All prefixes must be the same
        prefixes = set([bedline.amplicon_prefix for bedline in all_lines])
        if len(prefixes) != 1:
            raise ValueError(
                f"All bedlines must have the same prefix, ({','.join(prefixes)})"
            )
        self.prefix = prefixes.pop()

        # Check all chrom are the same
        chroms = set([bedline.chrom for bedline in all_lines])
        if len(chroms) != 1:
            raise ValueError(
                f"All bedlines must be on the same chromosome, ({','.join(chroms)})"
            )
        self.chrom = chroms.pop()
        # Check all pools are the same
        pools = set([bedline.pool for bedline in all_lines])
        if len(pools) != 1:
            raise ValueError(
                f"All bedlines must be in the same pool, ({','.join(map(str, pools))})"
            )
        self.pool = pools.pop()
        # Check all amplicon numbers are the same
        amplicon_numbers = set([bedline.amplicon_number for bedline in all_lines])
        if len(amplicon_numbers) != 1:
            raise ValueError(
                f"All bedlines must be the same amplicon, ({','.join(map(str, amplicon_numbers))})"
            )
        self.amplicon_number = amplicon_numbers.pop()

    @property
    def ipool(self) -> int:
        """Return the 0-based pool number"""
        return self.pool - 1

    def merge(self) -> tuple[BedLine, BedLine]:
        """
        Merges multiple forward and reverse primers into a single bedline
        """
        # Handle the fbedlines
        fbedline_start = min([bedline.start for bedline in self.fbedlines])
        fbedline_end = max([bedline.end for bedline in self.fbedlines])
        # Give it the longest seq to maintain the 7col format
        fbedline_sequence = max(
            [bedline.sequence for bedline in self.fbedlines], key=len
        )

        rbedline_start = min([bedline.start for bedline in self.rbedlines])
        rbedline_end = max([bedline.end for bedline in self.rbedlines])
        rbedline_sequence = max(
            [bedline.sequence for bedline in self.rbedlines], key=len
        )

        return BedLine(
            chrom=self.chrom,
            start=fbedline_start,
            end=fbedline_end,
            primername=f"{self.prefix}_{self.amplicon_number}_LEFT_1",
            pool=self.pool,
            strand=StrandEnum.FORWARD,
            sequence=fbedline_sequence,
        ), BedLine(
            chrom=self.chrom,
            start=rbedline_start,
            end=rbedline_end,
            primername=f"{self.prefix}_{self.amplicon_number}_RIGHT",
            pool=self.pool,
            strand=StrandEnum.REVERSE,
            sequence=rbedline_sequence,
        )
