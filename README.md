# PrimalBedTools

PrimalBedTools is a library for manipulating and processing BED files, particularly focused on primer-related operations. It provides several functions for common BED file operations including coordinate remapping, sorting, updating, and amplicon generation.

Functions are wrapped in a CLI for ease of use.

## Installation

```bash
pip install primalbedtools
```

## Usage

```bash
primalbedtools <command> [options]
```

### Commands

#### Remap
Remaps coordinates in a BED file using a multiple sequence alignment (MSA).

```bash
primalbedtools remap --bed input.bed --msa alignment.fasta --from_id source_id --to_id target_id
```

Required arguments:
- `--bed`: Input BED file
- `--msa`: Multiple sequence alignment file in FASTA format
- `--from_id`: Source sequence ID to remap from
- `--to_id`: Target sequence ID to remap to

#### Sort
Sorts entries in a BED file.

```bash
primalbedtools sort input.bed
```

Required arguments:
- `bed`: Input BED file to sort

#### Update
Updates information in a BED file, particularly primer names.

```bash
primalbedtools update input.bed
```

Required arguments:
- `bed`: Input BED file to update

#### Amplicon
Creates an amplicon BED file from primer pairs.

```bash
primalbedtools amplicon input.bed [--primertrim]
```

Required arguments:
- `bed`: Input BED file containing primer information

Optional arguments:
- `-t, --primertrim`: Generate primer-trimmed amplicons

## Examples

1. Remap coordinates from one reference to another:
```bash
primalbedtools remap --bed primers.bed --msa refs.fasta --from_id ref1 --to_id ref2
```

2. Sort a BED file:
```bash
primalbedtools sort input.bed > sorted.bed
```

3. Generate amplicons with primer trimming:
```bash
primalbedtools amplicon primers.bed --primertrim > trimmed_amplicons.bed
```

## Output

All commands output modified BED files to stdout while preserving the original header information. The amplicon command outputs either standard or primer-trimmed amplicon coordinates based on the provided options.
