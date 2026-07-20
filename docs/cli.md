# PrimalBedTools CLI

PrimalBedTools provides a command-line interface for working with primer BED files. It offers utilities for manipulating, validating, and analyzing primer designs.

## Installation

Install primalbedtools using pip:

```bash
pip install primalbedtools
```

or conda:

```bash
conda install bioconda::primalbedtools
```

## Basic Usage

```bash
primalbedtools <command> [options]
```

## Global options

These apply to every command and must be given **before** the command name.

- `-v`, `--verbose`: Enable debug-level logging
- `-q`, `--quiet`: Suppress warnings, leaving only errors

```bash
primalbedtools -v amplicon primers.bed   # correct
primalbedtools amplicon -v primers.bed   # error
```

Warnings and other advisory messages are written to stderr, so redirecting
stdout captures only the data:

```bash
primalbedtools amplicon primers.bed > amplicons.bed   # warnings still shown
primalbedtools amplicon primers.bed > amplicons.bed 2>/dev/null   # silenced
primalbedtools amplicon primers.bed > amplicons.bed 2>&1   # folded into the file
```

## Commands

### remap

Remap BED file coordinates from one reference to another using a multiple sequence alignment.

```bash
primalbedtools remap --bed <bed_file> --msa <msa_file> --from_id <source_id> --to_id <target_id>
```

**Arguments:**

- `--bed`: Input BED file (required)
- `--msa`: Multiple sequence alignment file (required)
- `--from_id`: Source sequence ID to remap from (required)
- `--to_id`: Target sequence ID to remap to (required)

**Example:**
```bash
primalbedtools remap --bed primers.bed --msa alignment.fasta --from_id MN908947.3 --to_id BA.2
```

### sort

Sort BED file by chromosome, amplicon number, and primer direction.

```bash
primalbedtools sort <bed_file>
```

**Arguments:**

- `bed`: Input BED file

**Example:**
```bash
primalbedtools sort primers.bed > primers.sorted.bed
```

### update

Update primer names to v2 format (`prefix_number_DIRECTION_index`).

```bash
primalbedtools update <bed_file>
```

**Arguments:**

- `bed`: Input BED file

**Example:**
```bash
primalbedtools update primers.v1.bed > primers.v2.bed
```

### amplicon

Generate amplicon information from primer pairs.

```bash
primalbedtools amplicon <bed_file> [--primertrim]
```

**Arguments:**

- `bed`: Input BED file
- `-t, --primertrim`: Generate primer-trimmed amplicon information

**Example:**
```bash
primalbedtools amplicon primers.bed > amplicons.txt
primalbedtools amplicon primers.bed --primertrim > trimmed_amplicons.txt
```

**Mixed prefixes:** primers belonging to one amplicon normally share a prefix.
Where they don't, the amplicon is named using every prefix joined with `-`, so
the mismatch stays visible in the output rather than being silently resolved to
one of them:

```
SARS-CoV-2-F_1_LEFT_1
SARS-CoV-2-R_1_RIGHT_1   ->   SARS-CoV-2-F-SARS-CoV-2-R_1
```

A single warning summarising every affected amplicon is written to stderr. Use
`-v` to list them individually, or `-q` to suppress it.

### merge

Merge primers with the same properties (chromosome, amplicon number, direction).

```bash
primalbedtools merge <bed_file>
```

**Arguments:**

- `bed`: Input BED file

**Example:**
```bash
primalbedtools merge primers.bed > primers.merged.bed
```

### fasta

Convert BED file to FASTA format.

```bash
primalbedtools fasta <bed_file>
```

**Arguments:**

- `bed`: Input BED file

**Example:**
```bash
primalbedtools fasta primers.bed > primers.fasta
```

### validate_bedfile

Validate a BED file for internal consistency (correct primer pairings, etc.).

```bash
primalbedtools validate_bedfile <bed_file>
```

**Arguments:**

- `bed`: Input BED file

**Example:**
```bash
primalbedtools validate_bedfile primers.bed
```

### validate

Validate a BED file against a reference genome.

```bash
primalbedtools validate <bed_file> <fasta_file>
```

**Arguments:**

- `bed`: Input BED file
- `fasta`: Reference FASTA file

**Example:**
```bash
primalbedtools validate primers.bed reference.fasta
```

### downgrade

Downgrade a BED file from v2 to v1 primer name format.

```bash
primalbedtools downgrade <bed_file> [--merge-alts]
```

**Arguments:**

- `bed`: Input BED file
- `--merge-alts`: Merge alternative primers (removes _alt suffixes)

**Example:**
```bash
# Downgrade with alternative primers
primalbedtools downgrade primers.v2.bed > primers.v1.bed

# Downgrade without alternative primers
primalbedtools downgrade primers.v2.bed --merge-alts > primers.v1.merged.bed
```
