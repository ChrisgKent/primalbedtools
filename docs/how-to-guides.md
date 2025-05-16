## Reading in bed files

The first step in using bedfiles in reading them into the BedLine class. This can be done from either a file or from a string.

```python
from primalbedtools.bedfiles import BedLineParser, BedLine

# From a file
header, bedlines = BedLineParser.from_file('./primer.bed')

# From a str
header, bedlines = BedLineParser.from_str('# header!!!\nMN908947.3\t47\t78\tSARS-CoV-2_1_LEFT_1\t1\t+\tCTCTTGTAGATCTGTTCTCTAAACGAACTTT\nMN908947.3\t419\t447\tSARS-CoV-2_1_RIGHT_1\t1\t-\tAAAACGCCTTTTTCAACTTCTACTAAGC\n')

# Headers 
print(header)
['# header!!!']

# Bedlines 
print([bl.primername for bl in bedlines])
['SARS-CoV-2_1_LEFT_1', 'SARS-CoV-2_1_RIGHT_1']

```

## Writing bedfiles

BedLines can be converted into a string or written to a file
```python
# Create the str
bed_str = BedLineParser.to_str(headers = ['example header'], bedlines = bedlines)

# Write to a file
BedLineParser.to_file("./primer.bed", headers = [], bedlines = bedlines)
```

## Using the BedLine

The BedLine class holds all the information contained in the bed file, while also providing field validation (ensuring the BedLine is valid at all times).

### Properties 

All expected properties are accessible.  Including some calculated ones. Please see [BedLine reference](reference.md)

```python
bl = bedlines[0]

print(bl.start)
47

print(bl.primername)
'SARS-CoV-2_1_LEFT_1'
```

### Setting Properties

Some properties are reliant on others. For example `primername` depends on `amplicon_prefix, amplicon_number, strand, primer_suffix`. Hence updating one will update all connected values. 

```python
# see old name
print(bl.primername)
'SARS-CoV-2_1_LEFT_1'

# update amplicon_prefix
bl.amplicon_prefix = "new-amp-name"

# see updated primername
print(bl.primername)
'new-amp-name_1_LEFT_1'

# change the strand
bl.strand = "-"
print(bl.primername)
'new-amp-name_1_RIGHT_1'

# this also works in reverse 
bl.primername = "final_10_LEFT_alt1"
print(bl.amplicon_prefix, bl.amplicon_number, bl.strand, bl.primer_suffix, sep= " | " )
"final | 10 | + | alt1"
```


### Internal validation 

Attempting to set an invalid property will raise an error.

```python
bl.primername = "test"
ValueError: Invalid primername: (test). Must be in v1 or v2 format

bl.pool = 0
ValueError: pool is 1-based pos int pool number. Got (0)

bl.start = "A"
ValueError: start must be an int. Got (A)

bl.amplicon_prefix = "test_invalid"
ValueError: Invalid amplicon_prefix: (test_invalid). Must be alphanumeric or hyphen.
```


## Common Operations

### Converting between versions

During development of primer.bed files multiple versions have been used. 

The current version (v2) uses a primername in the form of `{amplicon_prefix}_{amplicon_number}_{strand}_{primer_number}` (eg `SARS-CoV-2_1_RIGHT_1`). 

The old form typically used `{amplicon_prefix}_{amplicon_number}` or `{amplicon_prefix}_{amplicon_number}_alt1` (eg `SARS-CoV-2_1_RIGHT` | `SARS-CoV-2_1_RIGHT_alt`) to represent multiple primers for a site. 

Therefore being able to convert between the two is useful as certain tool product / expect different versions. 

#### Updating 

v1 bedfile example
```
MN908947.3	47	78	SARS-CoV-2_1_LEFT	1	+	CTCTTGTAGATCTGTTCTCTAAACGAACTTT
MN908947.3	419	447	SARS-CoV-2_1_RIGHT	1	-	AAAACGCCTTTTTCAACTTCTACTAAGC
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_alt	1	-	AAAACGCCTTTTTCAACTTCTACTAAAA
MN908947.3	344	366	SARS-CoV-2_2_LEFT	2	+	TCGTACGTGGCTTTGGAGACTC
MN908947.3	707	732	SARS-CoV-2_2_RIGHT	2	-	TCTTCATAAGGATCAGTGCCAAGCT
...
```

```python
from primalbedtools.bedfiles import BedFileModifier
updated_bl = BedFileModifier.update_primernames(bedlines)
```

updated to v2 
```
MN908947.3	47	78	SARS-CoV-2_1_LEFT_1	1	+	CTCTTGTAGATCTGTTCTCTAAACGAACTTT
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_2	1	-	AAAACGCCTTTTTCAACTTCTACTAAGC
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_1	1	-	AAAACGCCTTTTTCAACTTCTACTAAAA
MN908947.3	344	366	SARS-CoV-2_2_LEFT_1	2	+	TCGTACGTGGCTTTGGAGACTC
MN908947.3	707	732	SARS-CoV-2_2_RIGHT_1	2	-	TCTTCATAAGGATCAGTGCCAAGCT
...
```

#### Downgrading

!!! warning "The original order of the primers is lost. Therefore, combining upgrading and downgrading might produce different primernames (if unordered input)."

This will take the upgraded v2 output and convert it back into v1 
    
```python
from primalbedtools.bedfiles import BedFileModifier
downgraded_bl = BedFileModifier.downgrade_primernames(updated_bl)
```
downgraded to v1 (with alts)
```
MN908947.3	47	78	SARS-CoV-2_1_LEFT	1	+	CTCTTGTAGATCTGTTCTCTAAACGAACTTT
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_alt1 1	-	AAAACGCCTTTTTCAACTTCTACTAAGC
MN908947.3	419	447	SARS-CoV-2_1_RIGHT	1	-	AAAACGCCTTTTTCAACTTCTACTAAAA
MN908947.3	344	366	SARS-CoV-2_2_LEFT	2	+	TCGTACGTGGCTTTGGAGACTC
MN908947.3	707	732	SARS-CoV-2_2_RIGHT	2	-	TCTTCATAAGGATCAGTGCCAAGCT
...
```
!!! note "In this example the `alt` suffix changed to `alt1` moved from the primer with the sequence `AAAACGCCTTTTTCAACTTCTACTAAAA` to the primer with sequence `AAAACGCCTTTTTCAACTTCTACTAAGC`"

The presence of `_alt` can still provide breaking changes for some tools, and can be removed with `merge_alts`.

```python
from primalbedtools.bedfiles import BedFileModifier
downgraded_bl = BedFileModifier.downgrade_primernames(updated_bl, merge_alts=True)
```
downgraded to v1 (no alts)
```
MN908947.3	47	78	SARS-CoV-2_1_LEFT	1	+	CTCTTGTAGATCTGTTCTCTAAACGAACTTT
MN908947.3	419	447	SARS-CoV-2_1_RIGHT	1	-	AAAACGCCTTTTTCAACTTCTACTAAAA
MN908947.3	344	366	SARS-CoV-2_2_LEFT	2	+	TCGTACGTGGCTTTGGAGACTC
MN908947.3	707	732	SARS-CoV-2_2_RIGHT	2	-	TCTTCATAAGGATCAGTGCCAAGCT
...
```

### Sorting bedfiles

!!! note "Sorting a primer.bed file is having bedlines in amplicon order rather than index order."

```
MN908947.3	344	366	SARS-CoV-2_2_LEFT_1	2	+	TCGTACGTGGCTTTGGAGACTC
MN908947.3	47	78	SARS-CoV-2_1_LEFT_1	1	+	CTCTTGTAGATCTGTTCTCTAAACGAACTTT
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_2	1	-	AAAACGCCTTTTTCAACTTCTACTAAGC
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_1	1	-	AAAACGCCTTTTTCAACTTCTACTAAAA
MN908947.3	707	732	SARS-CoV-2_2_RIGHT_1	2	-	TCTTCATAAGGATCAGTGCCAAGCT
```

This will order by `chrom`, `amplicon_number`, `strand`, and `primer_suffix`.

```python
from primalbedtools.bedfiles import BedFileModifier
sorted_bl = BedFileModifier.sort_bedlines(bedlines)
```
sorted bedlines
```
MN908947.3	47	78	SARS-CoV-2_1_LEFT_1	1	+	CTCTTGTAGATCTGTTCTCTAAACGAACTTT
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_2	1	-	AAAACGCCTTTTTCAACTTCTACTAAGC
MN908947.3	419	447	SARS-CoV-2_1_RIGHT_1	1	-	AAAACGCCTTTTTCAACTTCTACTAAAA
MN908947.3	344	366	SARS-CoV-2_2_LEFT_1	2	+	TCGTACGTGGCTTTGGAGACTC
MN908947.3	707	732	SARS-CoV-2_2_RIGHT_1	2	-	TCTTCATAAGGATCAGTGCCAAGCT
```

### Grouping Bedlines

A common operation is to group bedlines by a property. This can be done with `group_by_chrom`, `group_by_amplicon_number`, `group_by_strand`, `group_primer_pairs`. 

```python
from primalbedtools.bedfiles import group_by_chrom, group_by_amplicon_number, group_by_strand, group_primer_pairs, group_by_pool

group_by_chrom(bedlines)
{'MN908947.3': [...]}

group_by_amplicon_number(bedlines)
{1: [...],
 2: [...]}

group_by_strand(bedlines)
{'+': [...],
 '-': [...]}

 group_by_pool(bedlines)
{1: [...],
 2: [...]}

# For each amplicon. It will return ([f bedlines], [r bedlines])
group_primer_pairs(bedlines)
[([<primalbedtools.bedfiles.BedLine object at 0x105431650>], [<primalbedtools.bedfiles.BedLine object at 0x105431fd0>, <primalbedtools.bedfiles.BedLine object at 0x105432850>]), ([<primalbedtools.bedfiles.BedLine object at 0x105432150>], [<primalbedtools.bedfiles.BedLine object at 0x105432650>])]
```