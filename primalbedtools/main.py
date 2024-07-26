import enum

from pydantic import BaseModel


class StrandEnum(enum.Enum):
    FORWARD = "+"
    NEGATIVE = "-"


class BedLine(BaseModel):
    """ """

    class Config:
        use_enum_values = True

    # properties
    chrom: str
    start: int
    end: int
    primername: str
    pool: int
    strand: StrandEnum
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

    def to_bed(self) -> str:
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.primername}\t{self.pool}\t{self.strand}\t{self.sequence}\n"


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


def read_bedfile(bedfile: str) -> tuple[list[str], list[BedLine]]:
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


def write_bedfile(bedfile: str, headers: list[str] | None, bedlines: list[BedLine]):
    with open(bedfile, "w") as f:
        f.write(create_bedfile_str(headers, bedlines))
