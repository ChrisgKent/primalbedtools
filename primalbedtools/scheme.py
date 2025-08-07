from typing import Optional

from primalbedtools.bedfiles import (
    BedLine,
    PrimerClass,
    bedline_from_str,
    create_bedfile_str,
    parse_headers_to_dict,
    read_bedfile,
    sort_bedlines,
    write_bedfile,
)


class Scheme:
    headers: list[str]
    bedlines: list[BedLine]

    def __init__(self, headers: Optional[list[str]], bedlines: list[BedLine]):
        # Parse the headers
        if headers is None:
            headers = []

        self.headers = headers
        self.bedlines = bedlines

    # io
    @classmethod
    def from_str(cls, str: str):
        headers, bedlines = bedline_from_str(str)
        return cls(headers, bedlines)

    @classmethod
    def from_file(cls, file: str):
        headers, bedlines = read_bedfile(file)
        return cls(headers, bedlines)

    def to_str(self) -> str:
        """returns the scheme as the bedfile string"""
        return create_bedfile_str(self.headers, self.bedlines)

    def to_file(self, path: str):
        """writes the bedline to the file path"""
        return write_bedfile(path, self.headers, self.bedlines)

    # modifiers
    def sort_bedlines(self):
        self.bedlines = sort_bedlines(self.bedlines)

    # properties
    @property
    def contains_probes(self) -> bool:
        """
        returns True if any of the bedlines are PROBES
        """
        for bedline in self.bedlines:
            if bedline.primer_class == PrimerClass.PROBE:
                return True
        return False

    @property
    def header_dict(self):
        return parse_headers_to_dict(self.headers)
