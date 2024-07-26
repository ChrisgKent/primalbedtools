import enum
import pathlib

from pydantic import BaseModel, ConfigDict


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
    @staticmethod
    def from_file(bedfile: str | pathlib.Path) -> tuple[list[str], list[BedLine]]:
        return read_bedfile(bedfile=bedfile)

    @staticmethod
    def from_str(bedfile_str: str) -> tuple[list[str], list[BedLine]]:
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
        return create_bedfile_str(headers, bedlines)

    @staticmethod
    def to_file(
        bedfile: str | pathlib.Path, headers: list[str] | None, bedlines: list[BedLine]
    ):
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
