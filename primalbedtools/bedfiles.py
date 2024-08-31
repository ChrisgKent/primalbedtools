import enum
import pathlib
import re

from pydantic import BaseModel, ConfigDict, field_validator

# Regular expressions for primer names
V1_PRIMERNAME = r"^[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)(_ALT[0-9]*|_alt[0-9]*)*$"
V2_PRIMERNAME = r"^[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)_[0-9]+$"


class PrimerNameVersion(enum.Enum):
    INVALID = 0
    V1 = 1
    V2 = 2


def version_primername(primername: str) -> PrimerNameVersion:
    """
    Check the version of a primername.
    """
    if re.match(V1_PRIMERNAME, primername):
        return PrimerNameVersion.V1
    elif re.match(V2_PRIMERNAME, primername):
        return PrimerNameVersion.V2
    else:
        return PrimerNameVersion.INVALID


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
    model_config = ConfigDict(
        use_enum_values=True, str_strip_whitespace=True, validate_assignment=True
    )

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
    def validate_primername(cls, v):
        if version_primername(v) == PrimerNameVersion.INVALID:
            raise ValueError(f"Invalid primername: ({v}). Must be in v1 or v2 format")
        return v

    @field_validator("pool")
    @classmethod
    def check_pool(cls, v):
        if v < 1:
            raise ValueError("Pool number must be greater than 0")
        return v

    # calculated properties
    @property
    def length(self):
        return self.end - self.start

    @property
    def primername_version(self) -> PrimerNameVersion:
        return version_primername(self.primername)

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


def group_primer_pairs(
    bedlines: list[BedLine],
) -> list[tuple[list[BedLine], list[BedLine]]]:
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
                    strand_to_bedlines.get(StrandEnum.FORWARD.value, []),
                    strand_to_bedlines.get(StrandEnum.REVERSE.value, []),
                )
            )

    return primer_pairs


def update_primernames(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Update primer names to v2 format in place.
    """
    # group the bedlines into primerpairs
    primer_pairs = group_primer_pairs(bedlines)

    # Update the primer names
    for fbedlines, rbedlines in primer_pairs:
        # Sort the bedlines by sequence
        fbedlines.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(fbedlines, start=1):
            bedline.primername = (
                f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_LEFT_{i}"
            )

        rbedlines.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(rbedlines, start=1):
            bedline.primername = (
                f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_RIGHT_{i}"
            )

    return bedlines


def downgrade_primernames(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Downgrades primer names to v1 format in place.
    """
    # group the bedlines into primerpairs
    primer_pairs = group_primer_pairs(bedlines)

    # Update the primer names
    for fbedlines, rbedlines in primer_pairs:
        # Sort the bedlines by sequence
        fbedlines.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(fbedlines, start=1):
            alt = "" if i == 1 else f"_alt{i-1}"
            bedline.primername = (
                f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_LEFT{alt}"
            )

        rbedlines.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(rbedlines, start=1):
            alt = "" if i == 1 else f"_alt{i-1}"
            bedline.primername = (
                f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_RIGHT{alt}"
            )

    return bedlines


def sort_bedlines(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Sorts bedlines by chrom, start, end, primername.
    """
    primerpairs = group_primer_pairs(bedlines)
    primerpairs.sort(key=lambda x: (x[0][0].chrom, x[0][0].amplicon_number))
    return [
        bedline
        for fbedlines, rbedlines in primerpairs
        for bedline in fbedlines + rbedlines
    ]


class BedFileModifier:
    """
    Collection of methods for modifying BED files.
    """

    @staticmethod
    def update_primernames(
        bedlines: list[BedLine],
    ) -> list[BedLine]:
        """
        Update primer names to v2 format in place.
        """
        return update_primernames(bedlines)

    @staticmethod
    def sort_bedlines(
        bedlines: list[BedLine],
    ) -> list[BedLine]:
        """
        Sorts the bedlines by chrom, amplicon number, direction, and sequence.
        """
        return sort_bedlines(bedlines)
