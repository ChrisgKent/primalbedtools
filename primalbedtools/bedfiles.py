import enum
import pathlib
import re
import typing
from typing import Optional, Union

from primalbedtools.utils import expand_ambiguous_bases, rc_seq, strip_all_white_space

# Regular expressions for primer names
V1_PRIMERNAME = r"^[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT|PROBE)(_ALT[0-9]*|_alt[0-9]*)*$"
V2_PRIMERNAME = r"^[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT|PROBE)_[0-9]+$"

AMPLICON_PREFIX = r"^[a-zA-Z0-9\-]+$"  # any alphanumeric or hyphen

V1_PRIMER_SUFFIX = r"^(ALT[0-9]*|alt[0-9]*)*$"

CHROM_REGEX = r"^[a-zA-Z0-9_.]+$"

PRIMER_WEIGHT_KEY = "pw"


class PrimerNameVersion(enum.Enum):
    INVALID = 0
    V1 = 1
    V2 = 2


class PrimerClass(enum.Enum):
    LEFT = "LEFT"
    RIGHT = "RIGHT"
    PROBE = "PROBE"


class Strand(enum.Enum):
    FORWARD = "+"
    REVERSE = "-"


def primer_class_str_to_enum(s: str):
    s = s.upper().strip()
    if s == "LEFT":
        return PrimerClass.LEFT
    elif s == "RIGHT":
        return PrimerClass.RIGHT
    elif s == "PROBE":
        return PrimerClass.PROBE

    raise ValueError(
        f"unknown primer class str ({s}). Should be ['RIGHT', 'LEFT', 'PROBE']"
    )


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


def check_amplicon_prefix(amplicon_prefix: str) -> bool:
    """
    Check if an amplicon prefix is valid.
    """
    return bool(re.match(AMPLICON_PREFIX, amplicon_prefix))


def check_valid_class_and_strand(cls: str, strand: str) -> bool:
    """
    This takes the PrimerClass str (LEFT|RIGHT|PROBE) and the strand str (+|-)
    and returns True is compatible.
    """
    # Probe can be any
    if cls == PrimerClass.PROBE.value:
        return True
    # LEFT primers must be forwards
    elif cls == PrimerClass.LEFT.value and strand == Strand.FORWARD.value:
        return True
    # RIGHT primers must be Reverse
    elif cls == PrimerClass.RIGHT.value and strand == Strand.REVERSE.value:
        return True
    else:
        return False


def parse_headers_to_dict(headers: list[str]) -> dict[str, str]:
    """
    parses the header strings into a dict.
    - Removes the leading # and any padding white space.
    - splits the header line on the first '=', with the key, value being the left, or right|None
    """
    attr_dict = {}

    for header in headers:
        header = header.rstrip().lstrip()  # remove lr whitespace

        if header.startswith("#"):
            header = header[1:].lstrip()  # Remove # and any padding

        parts = header.split("=", 1)

        if len(parts) == 1:
            attr_dict[parts[0].rstrip()] = None
        else:
            attr_dict[parts[0].rstrip()] = parts[1]

    return attr_dict


def create_primername(
    amplicon_prefix: str,
    amplicon_number: int,
    primer_class: PrimerClass,
    primer_suffix: Union[str, int, None],
):
    """
    Creates an unvalidated primername string.
    """
    values = [amplicon_prefix, amplicon_number, primer_class.value, primer_suffix]
    return "_".join([str(x) for x in values if x is not None])


def parse_primer_attributes_str(v: str) -> Optional[dict[str, str]]:
    new_dict = {}

    white_space_removed = strip_all_white_space(v)
    if not white_space_removed:
        return None

    # Check if old style weight is being provided.
    try:
        float(white_space_removed)
        # Parse old weight into a dict
        white_space_removed = f"{PRIMER_WEIGHT_KEY}={white_space_removed}"
    except ValueError:
        pass

    white_space_removed_kv = white_space_removed.split(";")
    if len(white_space_removed_kv) < 1:
        return None

    for kv in white_space_removed_kv:
        # Skip empty lines
        if not kv:
            continue

        try:
            k, v = kv.split("=")
        except ValueError:
            raise ValueError(
                f"Invalid PrimerAttributes: ({kv}). Must be in form k=v"
            ) from None
        # Check empty values
        if not k or not v:
            raise ValueError(f"Malformed k=v pair: ({kv})")
        # check dupe
        if k in new_dict:
            raise ValueError(f"Duplicate PrimerAttributes Key: ({k})")

        new_dict[k] = v

    # Only return a dict with values
    return new_dict if new_dict else None


def create_primer_attributes_str(
    primer_attributes: Union[dict[str, Union[str, float]], dict[str, str], None],
) -> Optional[str]:
    """
    Parses the dict into the ';' separated str. Strips all whitespace

    Returns None on empty dict or None
    """
    if primer_attributes is None or not primer_attributes:
        return None
    return ";".join(
        f"{strip_all_white_space(k)}={strip_all_white_space(str(v))}"
        for k, v in primer_attributes.items()
    )


def lr_string_to_strand_char(s: str) -> str:
    """
    Convert a LEFT/RIGHT string to a Strand.
    """
    parsed_strand = s.upper().strip()

    if parsed_strand == "LEFT":
        return Strand.FORWARD.value
    elif parsed_strand == "RIGHT":
        return Strand.REVERSE.value
    else:
        raise ValueError(f"Invalid strand: {s}. Must be LEFT or RIGHT")


def strand_char_to_primer_class_str(s: str) -> str:
    if s == Strand.FORWARD.value:
        return PrimerClass.LEFT.value
    elif s == Strand.REVERSE.value:
        return PrimerClass.RIGHT.value
    else:
        raise ValueError(f"unknown strand char ({s})")


def validate_strand(strand: str) -> str:
    """
    Validates and returns the strand. Raises ValueError
    """
    try:
        v = str(strand).strip()
    except ValueError as e:
        raise ValueError(f"strand must be a str. Got ({strand})") from e

    if v not in {x.value for x in Strand}:
        raise ValueError(
            f"strand must be a str of ({[x.value for x in Strand]}). Got ({v})"
        )
    return v


def validate_amplicon_prefix(amplicon_prefix: str):
    """
    Validates and returns the amplicon_prefix. Raises ValueError
    """
    try:
        v = str(amplicon_prefix).strip()
    except ValueError as e:
        raise ValueError(
            f"amplicon_prefix must be a str. Got ({amplicon_prefix})"
        ) from e

    if check_amplicon_prefix(v):
        return v
    else:
        raise ValueError(
            f"Invalid amplicon_prefix: ({v}). Must be alphanumeric or hyphen."
        )


def validate_primer_suffix(
    primer_suffix: Union[str, int, None],
) -> Union[str, int, None]:
    if primer_suffix is None:
        return None
    try:
        # Check for int. (handles _0 format)
        int_v = int(primer_suffix)
    except (ValueError, TypeError) as e:
        # If int() fails _alt format expected
        if isinstance(primer_suffix, str):
            if not re.match(V1_PRIMER_SUFFIX, primer_suffix):
                raise ValueError(
                    f"Invalid V1 primer_suffix: ({primer_suffix}). Must be `alt[0-9]*` or `ALT[0-9]*`"
                ) from e

            return primer_suffix

        else:
            raise ValueError(
                f"Invalid primer_suffix: ({primer_suffix}). Must be `alt[0-9]*` or `[0-9]`"
            ) from e

    if int_v < 0:
        raise ValueError(
            f"primer_suffix must be greater than or equal to 0. Got ({primer_suffix})"
        )
    else:
        return int_v


def validate_primer_name(primername: str) -> tuple[str, str, str, Union[str, None]]:
    """
    Validates the structure of the primer Name, and returns the unvalidated components.
    (Amplicon_prefix, Amplicon_number, Primer_class, Primer_suffix)
    """

    primername = strip_all_white_space(primername)
    if version_primername(primername) == PrimerNameVersion.INVALID:
        raise ValueError(
            f"Invalid primername: ({primername}). Must be in v1 or v2 format"
        )

    parts = primername.split("_")
    if len(parts) == 3:
        parts.append(None)  # type: ignore

    return (parts[0], parts[1], parts[2], parts[3])


class BedLine:
    """A class representing a single line in a primer.bed file.

    BedLine stores and validates all attributes of a primer entry in a BED file,
    with support for primername version handling, position validation, and
    output formatting. It maintains internal consistency between related fields
    like strand and primer direction, and automatically parses complex primername
    formats.

    Attributes:
        chrom (str): Chromosome name, must match pattern [a-zA-Z0-9_.]+
        start (int): 0-based start position of the primer
        end (int): End position of the primer
        primername (str): Name of the primer in either format v1 or v2
        pool (int): 1-based pool number (use ipool for 0-based pool number)
        strand (str): Strand of the primer ("+" for forward, "-" for reverse)
        sequence (str): Sequence of the primer
        attributes (dict[str,str|float], None): Dict of primer attributes (e.g., primerWeights (pw) for rebalancing).

    Properties:
        length (int): Length of the primer (end - start)
        amplicon_number (int): Amplicon number extracted from primername
        amplicon_prefix (str): Amplicon prefix extracted from primername
        primer_suffix (int, str, None): Suffix of the primer (integer index or alt string)
        primername_version (PrimerNameVersion): Version of the primername format
        ipool (int): 0-based pool number (pool - 1)
        direction_str (str): Direction as string ("LEFT" or "RIGHT")

    Examples:
        >>> from primalbedtools.bedfiles import BedLine
        >>> bedline = BedLine(
        ...     chrom="chr1",
        ...     start=100,
        ...     end=120,
        ...     primername="scheme_1_LEFT_alt1",
        ...     pool=1,
        ...     strand="+",
        ...     sequence="ACGTACGTACGTACGTACGT",
        ... )
        >>> print(bedline.length)
        20
        >>> print(bedline.primername_version)
        PrimerNameVersion.V1
        >>> print(bedline.to_bed())
        chr1	100	120	scheme_1_LEFT_alt1	1	+	ACGTACGTACGTACGTACGT
        >>> bedline.amplicon_prefix = "new-scheme"
        >>> print(bedline.to_bed())
        chr1    100     120     new-scheme_1_LEFT_alt1  1       +       ACGTACGTACGTACGTACGT
    """

    __slots__ = (
        "_chrom",
        "_start",
        "_end",
        "_pool",
        "_strand",
        "_sequence",
        "_attributes",
        "_amplicon_prefix",
        "_amplicon_number",
        "_primer_class",
        "_primer_suffix",
    )
    # properties
    _chrom: str
    _start: int
    # primername is a calculated property
    _end: int
    _pool: int
    _strand: str
    _sequence: str

    # primerAttributes
    _attributes: dict[str, Union[str, float]]

    # primernames components
    _amplicon_prefix: str
    _amplicon_number: int
    _primer_class: PrimerClass
    _primer_suffix: Union[int, str, None]

    def __init__(
        self,
        chrom: str,
        start: Union[int, str],
        end: Union[int, str],
        primername: str,
        pool: Union[int, str],
        strand: str,
        sequence: str,
        attributes: Union[dict[str, Union[str, float]], str, None, float] = None,
    ) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.pool = pool
        self.sequence = sequence
        # Set attributes
        self.attributes = attributes  # Ensure the type matches the expected UnionType

        # use private fields to sidestep primer_class to strand validation
        self._strand = validate_strand(strand)
        (
            self.amplicon_prefix,
            self.amplicon_number,
            primer_class,
            self.primer_suffix,
        ) = validate_primer_name(primername)
        self._primer_class = primer_class_str_to_enum(primer_class)

        # Check the primer_class and strand are compatible
        if not check_valid_class_and_strand(self.primer_class_str, self.strand):
            raise ValueError(
                f"primername ({self.primername}) implies direction ({self.primer_class_str}), which is incompatible with ({strand})"
            )

    @property
    def chrom(self):
        """Return the chromosome of the primer"""
        return self._chrom

    @chrom.setter
    def chrom(self, v):
        v = "".join(str(v).split())  # strip all whitespace
        if re.match(CHROM_REGEX, v):
            self._chrom = v
        else:
            raise ValueError(f"chrom must match '{CHROM_REGEX}'. Got (`{v}`)")

    @property
    def start(self):
        """Return the start position of the primer"""
        return self._start

    @start.setter
    def start(self, v):
        try:
            v = int(v)
        except (ValueError, TypeError) as e:
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
        except (ValueError, TypeError) as e:
            raise ValueError(f"end must be an int. Got ({v})") from e
        if v < 0:
            raise ValueError(f"end must be greater than or equal to 0. Got ({v})")
        self._end = v

    @property
    def amplicon_number(self) -> int:
        """Return the amplicon number of the primer"""
        return self._amplicon_number

    @property
    def amplicon_name(self) -> str:
        """Return the amplicon name of the primer"""
        return f"{self.amplicon_prefix}_{self.amplicon_number}"

    @amplicon_number.setter
    def amplicon_number(self, v):
        try:
            v = int(v)
        except ValueError as e:
            raise ValueError(f"amplicon_number must be an int. Got ({v})") from e
        if v < 0:
            raise ValueError(
                f"amplicon_number must be greater than or equal to 0. Got ({v})"
            )

        self._amplicon_number = v

    @property
    def amplicon_prefix(self) -> str:
        """Return the amplicon_prefix of the primer"""
        return self._amplicon_prefix

    @amplicon_prefix.setter
    def amplicon_prefix(self, v):
        try:
            v = str(v).strip()
        except ValueError as e:
            raise ValueError(f"amplicon_prefix must be a str. Got ({v})") from e

        if check_amplicon_prefix(v):
            self._amplicon_prefix = v
        else:
            raise ValueError(
                f"Invalid amplicon_prefix: ({v}). Must be alphanumeric or hyphen."
            )

    @property
    def primer_suffix(self) -> Union[int, str, None]:
        """Return the primer_suffix of the primer"""
        return self._primer_suffix

    @primer_suffix.setter
    def primer_suffix(self, v: Union[int, str, None]):
        self._primer_suffix = validate_primer_suffix(v)

    @property
    def primer_class(self) -> PrimerClass:
        """Return the PrimerDirection enum"""
        return self._primer_class

    @property
    def primer_class_str(self) -> str:
        """Return the string representation of PrimerDirection"""
        return self._primer_class.value

    @primer_class.setter
    def primer_class(self, v: str):
        new_primer_class = primer_class_str_to_enum(v)
        if check_valid_class_and_strand(new_primer_class.value, self.strand):
            self._primer_class = new_primer_class
        else:
            raise ValueError(
                f"The new primer_class ({new_primer_class.value}) is incompatible with current strand ({self.strand}). Please use method 'force_change' to update both."
            )

    @property
    def primername(self):
        """Return the primername of the primer"""
        return create_primername(
            self.amplicon_prefix,
            self.amplicon_number,
            self.primer_class,
            self.primer_suffix,
        )

    @primername.setter
    def primername(self, v, new_strand: Optional[str] = None):
        v = v.strip()
        if version_primername(v) == PrimerNameVersion.INVALID:
            raise ValueError(f"Invalid primername: ({v}). Must be in v1 or v2 format")

        # Parse the primername
        parts = v.split("_")
        self.amplicon_prefix = parts[0]
        self.amplicon_number = int(parts[1])

        # Get primer_class
        new_primer_class = primer_class_str_to_enum(parts[2])
        implied_strand: Optional[str] = None

        # if non-ambiguous set strand
        if new_primer_class != PrimerClass.PROBE:
            if new_primer_class == PrimerClass.LEFT:
                implied_strand = Strand.FORWARD.value
            elif new_primer_class == PrimerClass.RIGHT:
                implied_strand = Strand.REVERSE.value
            else:
                raise ValueError("Unknown strand!")

        # Set the strand and class
        if new_strand:
            self.force_change(new_primer_class.value, new_strand)
        elif implied_strand:
            self.force_change(new_primer_class.value, implied_strand)

        # Try to parse the primer_suffix
        try:
            self.primer_suffix = parts[3]
        except IndexError:
            self.primer_suffix = None

    @property
    def pool(self):
        """Return the 1-based pool number of the primer"""
        return self._pool

    @pool.setter
    def pool(self, v):
        try:
            v = int(v)
        except (ValueError, TypeError) as e:
            raise ValueError(f"pool must be an int. Got ({v})") from e
        if v < 1:
            raise ValueError(f"pool is 1-based pos int pool number. Got ({v})")
        self._pool = v

    @property
    def strand(self):
        """Return the strand of the primer"""
        return self._strand

    @property
    def strand_class(self) -> Strand:
        """Return the strand class of the primer"""
        if self.strand == "+":
            return Strand.FORWARD
        elif self.strand == "-":
            return Strand.REVERSE
        else:
            raise ValueError(
                f"Unknown strand value ({self.strand}) in {self.primername}"
            )

    @strand.setter
    def strand(self, v):
        new_s = validate_strand(v)

        if check_valid_class_and_strand(self.primer_class_str, new_s):
            self._strand = new_s
        else:
            raise ValueError(
                f"The new stand ({new_s}) is incompatible with current primer_class ({self.primer_class_str}). Please use method 'force_change' to update both."
            )

    @property
    def sequence(self):
        """Return the sequence of the primer"""
        return self._sequence

    @sequence.setter
    def sequence(self, v):
        if not isinstance(v, str):
            raise ValueError(f"sequence must be a str. Got ({v})")
        self._sequence = v.upper()

    @property
    def attributes(self):
        """Returns the primer attributes"""
        return self._attributes

    @attributes.setter
    def attributes(
        self, v: Optional[Union[dict[str, Union[str, float]], str, float]]
    ) -> None:
        # Parse string
        if isinstance(v, str):
            new_dict = parse_primer_attributes_str(v)
        # parse dict
        elif isinstance(v, dict):
            new_dict = v
        elif v is None:
            self._attributes = {}
            return
        else:
            raise ValueError(f"Invalid primer attributes. Got ({v})")

        if new_dict is None:
            self._attributes = {}
            return

        # Parse the new dict
        parsed_dict: dict[str, Union[str, float]] = {
            strip_all_white_space(str(k)): strip_all_white_space(str(v))
            for k, v in new_dict.items()
        }

        self._attributes = parsed_dict

        # Call to primer weight setter to validate
        if PRIMER_WEIGHT_KEY in parsed_dict:
            self.weight = parsed_dict[PRIMER_WEIGHT_KEY]

    @property
    def weight(self):
        """Return the weight of the primer from the attributes dict"""
        if self._attributes is None:
            return None
        return self._attributes.get(PRIMER_WEIGHT_KEY)

    @weight.setter
    def weight(self, v: Union[float, str, int, None]):
        # Catch Empty and None
        if v is None or v == "":
            if self._attributes is None:
                return
            # remove pw from the attributes
            self._attributes.pop(PRIMER_WEIGHT_KEY, None)
            return
        try:
            v = float(v)
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"weight must be a float, None or empty str (''). Got ({v})"
            ) from e
        if v < 0:
            raise ValueError(f"weight must be greater than or equal to 0. Got ({v})")

        # Set primerWeights as pw in self._attributes
        if self._attributes is None:
            self._attributes = {}
        self._attributes[PRIMER_WEIGHT_KEY] = v

    # methods
    def force_change(self, primer_class: str, strand: str):
        """
        Updates compatible primerclass and strand simultaneously.
        """
        new_primer_class = primer_class_str_to_enum(primer_class)
        new_s = validate_strand(strand)

        if check_valid_class_and_strand(new_primer_class.value, new_s):
            self._primer_class = new_primer_class
            self._strand = new_s
        else:
            raise ValueError(
                f"The new primer_class ({new_primer_class.value}) is incompatible with current strand ({self.strand})."
            )

    # calculated properties
    @property
    def length(self):
        """Return the length of the primer"""
        return self.end - self.start

    @property
    def primername_version(self) -> PrimerNameVersion:
        return version_primername(self.primername)

    @property
    def ipool(self) -> int:
        """Return the 0-based pool number"""
        return self.pool - 1

    @property
    def direction_str(self) -> str:
        """Return 'LEFT' or 'RIGHT' based on strand"""
        return "LEFT" if self.strand == Strand.FORWARD.value else "RIGHT"

    def to_bed(self) -> str:
        """Convert the BedLine object to a BED formatted string."""
        # If a attributes is provided print. Else print empty string

        attribute_str = create_primer_attributes_str(self.attributes)
        if attribute_str is None:
            attribute_str = ""
        else:
            attribute_str = "\t" + attribute_str
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.primername}\t{self.pool}\t{self.strand}\t{self.sequence}{attribute_str}\n"

    def to_fasta(self, rc=False) -> str:
        """Convert the BedLine object to a FASTA formatted string."""
        if rc:
            return f">{self.primername}-rc\n{rc_seq(self.sequence)}\n"
        return f">{self.primername}\n{self.sequence}\n"


class BedLineParser:
    """Collection of methods for BED file input/output operations.

    This class provides static methods for parsing and writing BED files
    in various formats, handling both file system operations and string
    conversions.
    """

    @staticmethod
    def from_file(
        bedfile: typing.Union[str, pathlib.Path],
    ) -> tuple[list[str], list[BedLine]]:
        """Reads and parses a BED file from disk.

        Reads a BED file from the specified path and returns the header lines
        and parsed BedLine objects.

        Args:
            bedfile: Path to the BED file. Can be a string or Path object.

        Returns:
            tuple: A tuple containing: `list[str]`: Header lines from the BED file (lines starting with '#') `list[BedLine]`: Parsed BedLine objects from the file

        Raises:
            FileNotFoundError: If the specified file doesn't exist
            ValueError: If the file contains invalid BED entries

        Examples:
            >>> from primalbedtools.bedfiles import BedLineParser
            >>> headers, bedlines = BedLineParser.from_file("primers.bed")
            >>> print(f"Found {len(bedlines)} primer entries")
        """
        return read_bedfile(bedfile=bedfile)

    @staticmethod
    def from_str(bedfile_str: str) -> tuple[list[str], list[BedLine]]:
        """Parses a BED file from a string.

        Parses a string containing BED file content and returns the header lines
        and parsed BedLine objects.

        Args:
            bedfile_str: String containing BED file content.

        Returns:
            tuple: A tuple containing:
                - list[str]: Header lines from the BED string (lines starting with '#')
                - list[BedLine]: Parsed BedLine objects from the string

        Raises:
            ValueError: If the string contains invalid BED entries
        """
        return bedline_from_str(bedfile_str)

    @staticmethod
    def to_str(headers: typing.Optional[list[str]], bedlines: list[BedLine]) -> str:
        """Creates a BED file string from headers and BedLine objects.

        Combines header lines and BedLine objects into a properly formatted
        BED file string.

        Args:
            headers: List of header strings (with or without leading '#')
                    or None for no headers.
            bedlines: List of BedLine objects to format as strings.

        Returns:
            str: A formatted BED file string with headers and entries.

        Examples:
            >>> from primalbedtools.bedfiles import BedLine, BedLineParser
            >>> bedlines = [BedLine(...)]  # List of BedLine objects
            >>> headers = ["Track name=primers"]
            >>> bed_string = BedLineParser.to_str(headers, bedlines)
        """
        return create_bedfile_str(headers, bedlines)

    @staticmethod
    def to_file(
        bedfile: typing.Union[str, pathlib.Path],
        headers: typing.Optional[list[str]],
        bedlines: list[BedLine],
    ) -> None:
        """Writes headers and BedLine objects to a BED file.

        Creates or overwrites a BED file at the specified path with
        the provided headers and BedLine objects.

        Args:
            bedfile: Path where the BED file will be written.
                    Can be a string or Path object.
            headers: List of header strings (with or without leading '#')
                    or None for no headers.
            bedlines: List of BedLine objects to write to the file.

        Raises:
            IOError: If the file cannot be written
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
        if len(bedline) < 8:
            attributes = None
        else:
            attributes = bedline[7]

        return BedLine(
            chrom=bedline[0],
            start=bedline[1],
            end=bedline[2],
            primername=bedline[3],
            pool=bedline[4],
            strand=bedline[5],
            sequence=bedline[6],
            attributes=attributes,
        )
    except IndexError as a:
        raise IndexError(
            f"Invalid BED line value: ({bedline}): has incorrect number of columns"
        ) from a


def bedline_from_str(bedline_str: str) -> tuple[list[str], list[BedLine]]:
    """
    Create a list of BedLine objects from a BED string.
    """
    headers = []
    bedlines = []
    for line in bedline_str.strip().split("\n"):
        line = line.strip()

        # Handle headers
        if line.startswith("#"):
            headers.append(line)
        elif line:
            bedlines.append(create_bedline(line.split("\t")))

    return headers, bedlines


def read_bedfile(
    bedfile: typing.Union[str, pathlib.Path],
) -> tuple[list[str], list[BedLine]]:
    with open(bedfile) as f:
        text = f.read()
        return bedline_from_str(text)


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
    """Groups a list of BedLine objects by chromosome.

    Takes a list of BedLine objects and organizes them into a dictionary
    where keys are chromosome names and values are lists of BedLine objects
    that belong to that chromosome.

    Args:
        list_bedlines: A list of BedLine objects to group.

    Returns:
        dict[str, list[BedLine]]: A dictionary mapping chromosome names (str)
            to lists of BedLine objects.

    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.chrom not in bedlines_dict:
            bedlines_dict[bedline.chrom] = []
        bedlines_dict[bedline.chrom].append(bedline)
    return bedlines_dict


def group_by_amplicon_number(list_bedlines: list[BedLine]) -> dict[int, list[BedLine]]:
    """Groups a list of BedLine objects by amplicon number.

    Takes a list of BedLine objects and organizes them into a dictionary
    where keys are amplicon numbers and values are lists of BedLine objects
    with that amplicon number.

    Args:
        list_bedlines: A list of BedLine objects to group.

    Returns:
        dict[int, list[BedLine]]: A dictionary mapping amplicon numbers (int)
            to lists of BedLine objects.

    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.amplicon_number not in bedlines_dict:
            bedlines_dict[bedline.amplicon_number] = []
        bedlines_dict[bedline.amplicon_number].append(bedline)
    return bedlines_dict


def group_by_strand(list_bedlines: list[BedLine]) -> dict[str, list[BedLine]]:
    """Groups a list of BedLine objects by strand.

    Takes a list of BedLine objects and organizes them into a dictionary
    where keys are strand values ("+" or "-") and values are lists of
    BedLine objects on that strand.

    Args:
        list_bedlines: A list of BedLine objects to group.

    Returns:
        dict[str, list[BedLine]]: A dictionary mapping strand values (str)
            to lists of BedLine objects.

    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.strand not in bedlines_dict:
            bedlines_dict[bedline.strand] = []
        bedlines_dict[bedline.strand].append(bedline)
    return bedlines_dict


def group_by_class(
    list_bedlines: list[BedLine],
) -> dict[str, list[BedLine]]:
    """Groups a list of BedLine objects by primer class.

    Takes a list of BedLine objects and organizes them into a dictionary
    where keys are class values and values are lists of
    BedLine objects in that class.

    Args:
        list_bedlines: A list of BedLine objects to group.

    Returns:
        dict[str, list[BedLine]]: A dictionary mapping class values (str)
            to lists of BedLine objects.

    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.primer_class_str not in bedlines_dict:
            bedlines_dict[bedline.primer_class_str] = []
        bedlines_dict[bedline.primer_class_str].append(bedline)

    return bedlines_dict


def group_by_pool(
    list_bedlines: list[BedLine],
) -> dict[int, list[BedLine]]:
    """Groups a list of BedLine objects by pool number.

    Takes a list of BedLine objects and organizes them into a dictionary
    where keys are pool numbers and values are lists of BedLine objects
    with that pool number.

    Args:
        list_bedlines: A list of BedLine objects to group.

    Returns:
        dict[int, list[BedLine]]: A dictionary mapping pool numbers (int)
            to lists of BedLine objects.

    """
    bedlines_dict = {}
    for bedline in list_bedlines:
        if bedline.pool not in bedlines_dict:
            bedlines_dict[bedline.pool] = []
        bedlines_dict[bedline.pool].append(bedline)
    return bedlines_dict


def group_primer_pairs(
    bedlines: list[BedLine],
) -> list[tuple[list[BedLine], list[BedLine]]]:
    """Groups Primer BedLine objects into primer pairs by chromosome and amplicon number.

    This function takes a list of BedLine objects and groups them based on chromosome
    and amplicon number. For each group, it creates a tuple containing the forward
    (LEFT) primers as the first element and the reverse (RIGHT) primers as the
    second element.

    Args:
        bedlines: A list of BedLine objects to group into primer pairs.

    Returns:
        A list of tuples, where each tuple contains: (First element: List of forward (LEFT) primers, Second element: List of reverse (RIGHT) primers)

    """
    primer_pairs = []

    # Group by chrom
    for chrom_bedlines in group_by_chrom(bedlines).values():
        # Group by amplicon number
        for amplicon_number_bedlines in group_by_amplicon_number(
            chrom_bedlines
        ).values():
            # Generate primer pairs
            strand_to_bedlines = group_by_class(amplicon_number_bedlines)
            primer_pairs.append(
                (
                    strand_to_bedlines.get(PrimerClass.LEFT.value, []),
                    strand_to_bedlines.get(PrimerClass.RIGHT.value, []),
                )
            )

    return primer_pairs


def group_amplicons(
    bedlines: list[BedLine],
) -> list[dict[str, list[BedLine]]]:
    """Groups all (including PROBES) BedLine objects into amplicons by chromosome and amplicon number.

    This function takes a list of BedLine objects and groups them based on chromosome
    and amplicon number. For each amplicon number, it creates a dict containing the LEFT,
    PROBE, and RIGHT as the keys pointing to a list of corresponding Bedlines.

    Args:
        bedlines: A list of BedLine objects to group into primer pairs.

    Returns:
        A list of dicts, with the key being the primer class string, and value a list of BedLines.

    """
    primer_pairs = []

    # Group by chrom
    for chrom_bedlines in group_by_chrom(bedlines).values():
        # Group by amplicon number
        for amplicon_number_bedlines in group_by_amplicon_number(
            chrom_bedlines
        ).values():
            # Generate primer pairs
            primer_pairs.append(group_by_class(amplicon_number_bedlines))

    return primer_pairs


def update_primernames(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Update primer names to v2 format in place.
    """
    # group the bedlines into Amplicons
    primer_pairs = group_amplicons(bedlines)

    # Update the primer names
    for dicts in primer_pairs:
        # left primer
        lp = dicts.get(PrimerClass.LEFT.value, [])
        lp.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(lp, start=1):
            bedline.primername = f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_{PrimerClass.LEFT.value}_{i}"

        # probes
        pp = dicts.get(PrimerClass.PROBE.value, [])
        pp.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(pp, start=1):
            bedline.primername = f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_{PrimerClass.PROBE.value}_{i}"

        # right primers
        rp = dicts.get(PrimerClass.RIGHT.value, [])
        rp.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(rp, start=1):
            bedline.primername = f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_{PrimerClass.RIGHT.value}_{i}"

    return bedlines


def downgrade_primernames(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Downgrades primer names to v1 format in place.
    """
    # group the bedlines into Amplicons
    primer_pairs = group_primer_pairs(bedlines)

    # Update the primer names
    for left, right in primer_pairs:
        # Sort the bedlines by sequence
        left.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(left, start=1):
            alt = "" if i == 1 else f"_alt{i - 1}"
            bedline.primername = (
                f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_LEFT{alt}"
            )

        right.sort(key=lambda x: x.sequence)
        for i, bedline in enumerate(right, start=1):
            alt = "" if i == 1 else f"_alt{i - 1}"
            bedline.primername = (
                f"{bedline.amplicon_prefix}_{bedline.amplicon_number}_RIGHT{alt}"
            )

    return bedlines


def sort_bedlines(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Sorts bedlines by chrom, start, end, primername.
    """
    amplicons = group_amplicons(bedlines)
    amplicons.sort(
        key=lambda x: (
            x[PrimerClass.LEFT.value][0].chrom,
            x[PrimerClass.LEFT.value][0].amplicon_number,
        )
    )  # Uses left primers

    # Sorted list
    sorted_list = []

    for dicts in amplicons:
        # Left primers
        lp = dicts.get(PrimerClass.LEFT.value, [])
        lp.sort(
            key=lambda x: x.primer_suffix if x.primer_suffix is not None else x.sequence
        )
        sorted_list.extend(lp)

        # Probes
        pp = dicts.get(PrimerClass.PROBE.value, [])
        pp.sort(
            key=lambda x: x.primer_suffix if x.primer_suffix is not None else x.sequence
        )
        sorted_list.extend(pp)

        # Right Primers
        rp = dicts.get(PrimerClass.RIGHT.value, [])
        rp.sort(
            key=lambda x: x.primer_suffix if x.primer_suffix is not None else x.sequence
        )
        sorted_list.extend(rp)

    return sorted_list


def merge_primers(bedlines: list[BedLine]) -> list[BedLine]:
    """
    merges bedlines with the same chrom, amplicon number and class.
    """
    merged_bedlines = []

    for pdict in group_amplicons(bedlines):
        # Merge forward primers
        left = pdict.get(PrimerClass.LEFT.value, [])
        if left:
            fbedline_start = min([bedline.start for bedline in left])
            fbedline_end = max([bedline.end for bedline in left])
            fbedline_sequence = max([bedline.sequence for bedline in left], key=len)
            merged_bedlines.append(
                BedLine(
                    chrom=left[0].chrom,
                    start=fbedline_start,
                    end=fbedline_end,
                    primername=f"{left[0].amplicon_prefix}_{left[0].amplicon_number}_{PrimerClass.LEFT.value}_1",
                    pool=left[0].pool,
                    strand=Strand.FORWARD.value,
                    sequence=fbedline_sequence,
                )
            )

        # Merge probes
        probes = pdict.get(PrimerClass.PROBE.value, [])
        if probes:
            raise ValueError(
                f"can't merge schemes containing probes ({[p.primername for p in probes]})"
            )

        # Merge reverse primers
        right = pdict.get(PrimerClass.RIGHT.value, [])
        if right:
            rbedline_start = min([bedline.start for bedline in right])
            rbedline_end = max([bedline.end for bedline in right])
            rbedline_sequence = max([bedline.sequence for bedline in right], key=len)

            merged_bedlines.append(
                BedLine(
                    chrom=right[0].chrom,
                    start=rbedline_start,
                    end=rbedline_end,
                    primername=f"{right[0].amplicon_prefix}_{right[0].amplicon_number}_{PrimerClass.RIGHT.value}_1",
                    pool=right[0].pool,
                    strand=Strand.REVERSE.value,
                    sequence=rbedline_sequence,
                )
            )
    return merged_bedlines


def expand_bedlines(bedlines: list[BedLine]) -> list[BedLine]:
    """
    Expands ambiguous bases in the primer sequences to all possible combinations.
    """
    expanded_bedlines = []

    for bedline in bedlines:
        for expand_seq in expand_ambiguous_bases(bedline.sequence):
            expanded_bedlines.append(
                # Create a bunch of bedlines with the same name
                BedLine(
                    chrom=bedline.chrom,
                    start=bedline.start,
                    end=bedline.end,
                    primername=bedline.primername,
                    pool=bedline.pool,
                    strand=bedline.strand,
                    sequence=expand_seq,
                    attributes=bedline.attributes,
                )
            )
    # update the bedfile names
    expanded_bedlines = update_primernames(expanded_bedlines)
    return expanded_bedlines


class BedFileModifier:
    """Collection of methods for modifying BED files.

    This class provides static methods for common BED file operations such as
    updating primer names, sorting BED lines, and merging BED lines with the
    same characteristics.
    """

    @staticmethod
    def update_primernames(
        bedlines: list[BedLine],
    ) -> list[BedLine]:
        """Updates primer names to v2 format in place.

        Converts all primer names in the provided BedLine objects to the v2 format
        (prefix_number_class_index). Groups primers by chromosome and amplicon
        number, then updates each group.

        Args:
            bedlines: A list of BedLine objects to update.

        Returns:
            list[BedLine]: The updated list of BedLine objects with v2 format names.

        Examples:
            >>> from primalbedtools.bedfiles import BedLine, BedFileModifier
            >>> bedlines = [BedLine(...)]  # List of BedLine objects
            >>> updated = BedFileModifier.update_primernames(bedlines)
        """
        return update_primernames(bedlines)

    @staticmethod
    def downgrade_primernames(
        bedlines: list[BedLine], merge_alts=False
    ) -> list[BedLine]:
        """Downgrades primer names to v1 format in place.

        Converts all primer names in the provided BedLine objects to the v1 format
        (prefix_number_class_ALT). Groups primers by chromosome and amplicon
        number, then updates each group.

        Args:
            bedlines: A list of BedLine objects to downgrade.

        Returns:
            list[BedLine]: The updated list of BedLine objects with v1 format names.
        """
        if merge_alts:
            # Merge the alt primers
            bedlines = merge_primers(bedlines)
        # Downgrade the primer names
        return downgrade_primernames(bedlines)

    @staticmethod
    def sort_bedlines(
        bedlines: list[BedLine],
    ) -> list[BedLine]:
        """Sorts the bedlines by chrom, amplicon number, class, and sequence.

        Groups BedLine objects into primer pairs, sorts those pairs by chromosome
        and amplicon number, then returns a flattened list of the sorted BedLine objects.

        Args:
            bedlines: A list of BedLine objects to sort.

        Returns:
            list[BedLine]: A new list containing the sorted original BedLine objects.

        Examples:
            >>> from primalbedtools.bedfiles import BedLine, BedFileModifier
            >>> bedlines = [BedLine(...)]  # List of BedLine objects
            >>> sorted_lines = BedFileModifier.sort_bedlines(bedlines)
        """
        return sort_bedlines(bedlines)

    @staticmethod
    def merge_primers(
        bedlines: list[BedLine],
    ) -> list[BedLine]:
        """Merges bedlines with the same chrom, amplicon number and class.

        Groups BedLine objects into primer pairs, then for each forward and reverse
        group, creates a merged BedLine with:
        - The earliest start position
        - The latest end position
        - The longest sequence
        - The amplicon prefix, number, and pool from the first BedLine

        Args:
            bedlines: A list of BedLine objects to merge.

        Returns:
            list[BedLine]: A new list containing new merged BedLine objects.

        Examples:
            >>> from primalbedtools.bedfiles import BedLine, BedFileModifier
            >>> bedlines = [BedLine(...)]  # List of BedLine objects
            >>> merged_lines = BedFileModifier.merge_primers(bedlines)
        """
        return merge_primers(bedlines)
