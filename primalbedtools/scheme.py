from typing import Optional

from primalbedtools.bedfiles import (
    BedLine,
    PrimerClass,
    bedline_from_str,
    create_bedfile_str,
    merge_primers,
    parse_headers_to_dict,
    read_bedfile,
    sort_bedlines,
    update_primernames,
    write_bedfile,
)

DEFAULT_CSV_HEADERS = [
    "chrom",
    "start",
    "end",
    "primername",
    "pool",
    "strand",
    "sequence",
    "amplicon_prefix",
    "amplicon_number",
    "primer_class_str",
    "primer_suffix",
]


class Scheme:
    """A class representing a primer scheme with headers and primer bed lines.

    A Scheme contains both the header lines (comments) and the primer bed lines
    that define a complete primer scheme for amplicon sequencing or qPCR.

    Please use Scheme.from_str() or Scheme.from_file() for creation.

    Attributes:
        headers (list[str]): List of header/comment lines from the bed file
        bedlines (list[BedLine]): List of BedLine objects representing primers
    """

    headers: list[str]
    bedlines: list[BedLine]

    def __init__(self, headers: Optional[list[str]], bedlines: list[BedLine]):
        """Initialize a Scheme with headers and bedlines.

        Please use Scheme.from_str() or Scheme.from_file() for creation.

        Args:
            headers: Optional list of header strings. If None, defaults to empty list.
            bedlines: List of BedLine objects representing the primers in the scheme.
        """
        # Parse the headers
        if headers is None:
            headers = []

        self.headers = headers
        self.bedlines = bedlines

    # io
    @classmethod
    def from_str(cls, str: str):
        """Create a Scheme from a bed file string.

        Args:
            str: String containing bed file content with headers and primer lines.

        Returns:
            Scheme: A new Scheme object parsed from the string.
        """
        headers, bedlines = bedline_from_str(str)
        return cls(headers, bedlines)

    @classmethod
    def from_file(cls, file: str):
        """Create a Scheme from a bed file on disk.

        Args:
            file: Path to the bed file to read.

        Returns:
            Scheme: A new Scheme object loaded from the file.
        """
        headers, bedlines = read_bedfile(file)
        return cls(headers, bedlines)

    def to_str(self) -> str:
        """Convert the scheme to a bed file format string.

        Returns:
            str: String representation of the scheme in bed file format,
                 including headers and all primer lines.
        """
        return create_bedfile_str(self.headers, self.bedlines)

    def to_file(self, path: str):
        """Write the scheme to a bed file on disk.

        Args:
            path: File path where the bed file should be written.
        """
        return write_bedfile(path, self.headers, self.bedlines)

    # modifiers
    def sort_bedlines(self):
        """Sort the bedlines in canonical order in place.

        Sorts bedlines by chromosome, amplicon number, direction, and primer suffix
        to ensure consistent ordering across the scheme.
        """
        self.bedlines = sort_bedlines(self.bedlines)

    def merge_primers(self):
        """merges bedlines with the same chrom, amplicon number and class in place"""
        self.bedlines = merge_primers(self.bedlines)

    def update_primernames(self):
        """Updates PrimerNames into V2 format in place"""
        self.bedlines = update_primernames(self.bedlines)

    # properties
    @property
    def contains_probes(self) -> bool:
        """Check if the scheme contains any probe primers.

        Returns:
            bool: True if any bedlines are PROBE type primers, False otherwise.
        """
        for bedline in self.bedlines:
            if bedline.primer_class == PrimerClass.PROBE:
                return True
        return False

    @property
    def header_dict(self) -> dict:
        """Parse headers into a dictionary format.

        Returns:
            dict: Dictionary representation of the header lines, parsed according
                  to common header formats used in bed files.
        """
        return parse_headers_to_dict(self.headers)

    def to_delim_str(
        self, include_headers: bool = True, use_header_aliases: bool = False
    ):
        return to_delim_str(
            self, include_headers=include_headers, use_header_aliases=use_header_aliases
        )


def to_delim_str(
    scheme: Scheme, include_headers: bool = True, use_header_aliases: bool = False
) -> str:
    """
    Turns a bedfile into a full expanded delim separated file
    """
    # Define the default headers
    headers = DEFAULT_CSV_HEADERS

    lines_to_write: list[str] = []

    header_aliases = scheme.header_dict
    aliases_to_attr = {v: k for k, v in header_aliases.items()}

    # Parse the attr strings add new headers
    for bl in scheme.bedlines:
        for k in bl.attributes.keys():
            if use_header_aliases:
                k = header_aliases.get(k, k)
            if k not in headers:
                headers.append(k)

    # Create a csv line for each bedline
    if include_headers:
        lines_to_write.append(",".join(headers))

    for bl in scheme.bedlines:
        bl_csv: list[str] = []
        for h in headers:
            r = None
            try:
                r = bl.__getattribute__(h)
            except AttributeError:
                # Search _attribute dict
                if h in bl.attributes:
                    r = bl.attributes[h]
                elif h in aliases_to_attr:
                    r = bl.attributes.get(aliases_to_attr[h])

            bl_csv.append(str(r) if r is not None else "")

        lines_to_write.append(",".join(bl_csv))

    # write all complete lines
    return "\n".join(lines_to_write)
