import enum
import pathlib
import re
import typing

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


class BedLine:
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

    # properties
    _chrom: str
    _start: int
    _end: int
    _primername: str
    _pool: int
    _strand: str
    _sequence: str

    def __init__(
        self,
        chrom,
        start,
        end,
        primername,
        pool,
        strand,
        sequence,
    ) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.primername = primername
        self.pool = pool
        self.strand = strand
        self.sequence = sequence

    @property
    def chrom(self):
        return self._chrom

    @chrom.setter
    def chrom(self, v):
        try:
            v = str(v)
        except ValueError as e:
            raise ValueError(f"chrom must be a str. Got ({v})") from e
        self._chrom = v

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, v):
        try:
            v = int(v)
        except ValueError as e:
            raise ValueError(f"start must be an int. Got ({v})") from e
        if v < 0:
            raise ValueError(f"start must be greater than or equal to 0. Got ({v})")
        self._start = v

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, v):
        try:
            v = int(v)
        except ValueError as e:
            raise ValueError(f"end must be an int. Got ({v})") from e
        if v < 0:
            raise ValueError(f"end must be greater than or equal to 0. Got ({v})")
        self._end = v

    @property
    def primername(self):
        return self._primername

    @primername.setter
    def primername(self, v):
        if version_primername(v) == PrimerNameVersion.INVALID:
            raise ValueError(f"Invalid primername: ({v}). Must be in v1 or v2 format")
        self._primername = v

    @property
    def pool(self):
        return self._pool

    @pool.setter
    def pool(self, v):
        try:
            v = int(v)
        except ValueError as e:
            raise ValueError(f"pool must be an int. Got ({v})") from e
        if v < 1:
            raise ValueError(f"pool is 1-based pos int pool number. Got ({v})")
        self._pool = v

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, v):
        try:
            v = str(v)
        except ValueError as e:
            raise ValueError(f"strand must be a str. Got ({v})") from e

        if v not in {x.value for x in StrandEnum}:
            raise ValueError(
                f"strand must be a str of ({[x.value for x in StrandEnum]}). Got ({v})"
            )
        self._strand = v

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, v):
        if not isinstance(v, str):
            raise ValueError(f"sequence must be a str. Got ({v})")
        self._sequence = v.upper()

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
    def from_file(
        bedfile: typing.Union[str, pathlib.Path],
    ) -> tuple[list[str], list[BedLine]]:
        """
        Read and parse a BED file and return a tuple of headers and BedLine objects.
        : param bedfile: typing.Union[str, pathlib.Path]
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
    def to_str(headers: typing.Optional[list[str]], bedlines: list[BedLine]) -> str:
        """
        Creates a BED string from the headers and BedLine objects.
        : param headers: typing.Optional[list[str]]
        : param bedlines: list[BedLine]
        : return: str
        """
        return create_bedfile_str(headers, bedlines)

    @staticmethod
    def to_file(
        bedfile: typing.Union[str, pathlib.Path],
        headers: typing.Optional[list[str]],
        bedlines: list[BedLine],
    ) -> None:
        """
        Creates a BED file from the headers and BedLine objects.
        : param bedfile: typing.Union[str, pathlib.Path]
        : param headers: typing.Optional[list[str]]
        : param bedlines: list[BedLine]
        """
        write_bedfile(bedfile, headers, bedlines)


def create_bedline(bedline: list[str]) -> BedLine:
    """
    Creates a BedLine object from a list of string values.

    :param bedline: list[str]
        A list of string values representing a BED line. The list should contain the following elements:
        - chrom: str, the chromosome name
        - start: str, the start position (will be converted to int)
        - end: str, the end position (will be converted to int)
        - primername: str, the name of the primer
        - pool: str, the pool number (will be converted to int)
        - strand: str, the strand ('+' or '-')
        - sequence: str, the sequence of the primer

    :return: BedLine
        A BedLine object created from the provided values.

    :raises ValueError:
        If any of the values cannot be converted to the appropriate type.
    :raises IndexError:
        If the provided list does not contain the correct number of elements.
    """
    try:
        return BedLine(
            chrom=bedline[0],
            start=bedline[1],
            end=bedline[2],
            primername=bedline[3],
            pool=bedline[4],
            strand=bedline[5],
            sequence=bedline[6],
        )
    except IndexError as a:
        raise IndexError(
            f"Invalid BED line value: ({bedline}): has incorrect number of columns"
        ) from a


def read_bedfile(
    bedfile: typing.Union[str, pathlib.Path],
) -> tuple[list[str], list[BedLine]]:
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


def create_bedfile_str(
    headers: typing.Optional[list[str]], bedlines: list[BedLine]
) -> str:
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
    bedfile: typing.Union[str, pathlib.Path],
    headers: typing.Optional[list[str]],
    bedlines: list[BedLine],
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
) -> dict[str, list[BedLine]]:
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


def merge_bedlines(bedlines: list[BedLine]) -> list[BedLine]:
    """
    merges bedlines with the same chrom, amplicon number and direction.
    """
    merged_bedlines = []

    for fbedlines, rbedlines in group_primer_pairs(bedlines):
        # Merge forward primers
        if fbedlines:
            fbedline_start = min([bedline.start for bedline in fbedlines])
            fbedline_end = max([bedline.end for bedline in fbedlines])
            fbedline_sequence = max(
                [bedline.sequence for bedline in fbedlines], key=len
            )
            merged_bedlines.append(
                BedLine(
                    chrom=fbedlines[0].chrom,
                    start=fbedline_start,
                    end=fbedline_end,
                    primername=f"{fbedlines[0].amplicon_prefix}_{fbedlines[0].amplicon_number}_LEFT_1",
                    pool=fbedlines[0].pool,
                    strand=StrandEnum.FORWARD.value,
                    sequence=fbedline_sequence,
                )
            )

        # Merge reverse primers
        if rbedlines:
            rbedline_start = min([bedline.start for bedline in rbedlines])
            rbedline_end = max([bedline.end for bedline in rbedlines])
            rbedline_sequence = max(
                [bedline.sequence for bedline in rbedlines], key=len
            )

            merged_bedlines.append(
                BedLine(
                    chrom=rbedlines[0].chrom,
                    start=rbedline_start,
                    end=rbedline_end,
                    primername=f"{rbedlines[0].amplicon_prefix}_{rbedlines[0].amplicon_number}_RIGHT_1",
                    pool=rbedlines[0].pool,
                    strand=StrandEnum.REVERSE.value,
                    sequence=rbedline_sequence,
                )
            )
    return merged_bedlines


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

    @staticmethod
    def merge_bedlines(
        bedlines: list[BedLine],
    ) -> list[BedLine]:
        """
        Merges bedlines with the same chrom, amplicon number and direction.
        """
        return merge_bedlines(bedlines)
