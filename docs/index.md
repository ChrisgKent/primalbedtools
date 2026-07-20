# Quick start

## Summary 
Primalbedtools is no dependency python library for the validation of primerschemes files (primer.bed + reference.fasta). It has been designed to carry out validation, common operation and file reading / writing. 


## History  
Over the years, the file that describes the primers for amplicon sequencing has undergone a lot of changes. However, they typically are text files, which each line representing a single primer in the PCR reaction.  

The initial version (informally **v1**) was utilised in 2016 to 2018, had 6 columns to describe the location of the primer but did not include the primer sequence. 

!!! warning "**v1** is now fully depreciated. We view primers sequences essential for reproducibility and therefore they are required"  

The next iteration (informally **v2**) added a 7th column which contained the primer sequence. It also provided some structure for __primerName__ (unique identifier for each primer in 4th column) `{schemeName|uuid}_{ampliconNumber}_{LEFT|RIGHT}` and an optional `{_alt}` to denote spike in primers. 

!!! info "**v2** has been superseded and should be updated to **v3**. For legacy uses `primalbedtools` can happily parse **v2** and even convert to **v3**"

The current generation (**v3.0.0**) has been formalised (see [here](https://github.com/artic-network/primerscheme-specs/blob/20816ff7cd53bdfaab7a605dee06de1c80be759f/pdf/primerscheme.pdf) for detailed specification) but to briefly summarise; 

- Probe based qPCR assays can be described with the `PROBE` class

- Comment lines starting with `#` are supported 

- An optional 8th col can include primer key-value metadata. For example a primer's gc content and score (ps) could be encoded as `gc=0.60;ps=100`

- primerNames must be in the form `{schemeName|uuid}_{ampliconNumber}_{LEFT|RIGHT|PROBE}_{primerNumber}`

!!! info "**v3** is the current file format. Most of ARTICnetwork's tools expect **v3** primer.bed files, but most use primalbedtools and hence should be able to use **v2**"

## Installation 

Install primalbedtools using pip:

```bash
pip install primalbedtools
```

or conda:

```bash
conda install bioconda::primalbedtools
```

or from source (requires uv)

```bash
git clone https://github.com/ChrisgKent/primalbedtools
cd primalbedtools 
uv sync 
uv run primalbedtools 
```

## Example bedfile

Here is a example bedfile, slightly modified for easy viewing. 
```
# chrom     start   end     primername              pool    strand  sequence            primerAttributes       
MN908947.3	47	    78	    SARS-CoV-2_1_LEFT_1	    1	    +	    CTCTTGTAGATCTT...   pw=1.0;ps=100
MN908947.3	419	    447	    SARS-CoV-2_1_RIGHT_1	1	    -	    AAAACGCCTTTCAA...   pw=0.8;ps=90
MN908947.3	344	    366	    SARS-CoV-2_2_LEFT_0	    2	    +	    TCGTACGTCTTTGG...   pw=1.0;ps=105
MN908947.3	707	    732	    SARS-CoV-2_2_RIGHT_0	2	    -	    TCTTCAAGGATCAG...   pw=1.2;ps=104
```

In primalbedtools each primer is represented by a BedLine object (see BedLine section for detailed docs).

The BedLine provides access to all expected fields.
```python
>>> bl.primername 
'SARS-CoV-2_1_LEFT_1'
>>> bl.sequence
'CTCTTGTAGATCTGTTCTCTAAACGAACTTT'
>>> bl.attributes
{'pw': 1.0, 'ps': '100'}
```

Alongside some calculated ones.
```python
>>> bl.length
31
>>> bl.primer_class_str
'LEFT'
>>> bl.amplicon_number
1
>>> bl.primer_suffix
1
>>> bl.ipool
0
```

