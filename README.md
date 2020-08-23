sequencetoolkit
===============
[![Build](https://github.com/informationsea/sequencetoolkit/workflows/Rust/badge.svg)](https://github.com/informationsea/sequencetoolkit/actions)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/informationsea/sequencetoolkit)](https://hub.docker.com/r/informationsea/sequencetoolkit)

Toolkit for Genome Sequence Analysis

Features
--------

### BED Utilities

1. Merge overlapped regions
2. Expand region size

### Gene Annotation Utilities

1. Convert CDS/Transcript position into genomic position
2. Convert genomic position into CDS/Transcript position

### VCF Utilities

1. Add `GenotypeCount` and `nhomalt` INFO tags.
2. Convert a VCF file to CSV, TSV and Excel file.
3. Replace contig names.
4. Replace sample names.
5. Remove INFO tags.
6. Remove FORMAT tags.
7. Remove non-standard headers.
8. Generate CREATE TABLE SQL from VCF file.

Usage
-----

### Add Genotype Count

Add or rewrite allele frequency, count, number and genotype count

```
USAGE:
    sequencetoolkit vcfutils add-af [OPTIONS] [--] [input]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --category <category>    category mapping file (csv or tsv)
    -i, --id <id>...             ID column name in category mapping file
    -o, --output <output>        Output CSV
    -v, --value <value>...       value column name in category mapping file

ARGS:
    <input>    Input VCF file
```

First line in a category mapping file is used as header. If you use category mapping file,
`id` and `value` options are required.

#### Example

* `category.csv`

```
id,sex,protocol
SAMPLE1,M,HiSeq2500
SAMPLE2,F,HiSeq2500
SAMPLE3,M,NovaSeq600
SAMPLE4,F,NovaSeq600
```

`vcfutils add-af --output output.vcf.gz --category category.csv --id id --value sex protocol -- input.vcf.gz`

Note: output vcf file is not compressed with bgzip.

### Convert a VCF file to CSV, TSV or Excel file

```
USAGE:
    sequencetoolkit vcfutils vcf2csv [FLAGS] [OPTIONS] [--] [input]

FLAGS:
    -d, --decode-genotype        Decode GT format tag number into alleles
    -h, --help                   Prints help information
    -m, --split-multi-allelic    Split multi allelic sites into multiple lines
    -V, --version                Prints version information

OPTIONS:
    -c, --canonical-list <canonical-list>    Canonical transcript list (created with extract-canonical command)
    -t, --data-type <datatype>               Data type to output [default: auto]  [possible values: auto,
                                             xlsx, csv, tsv]
    -f, --format <format>...                 FORMAT tags to include
    -i, --info <info>...                     INFO tags to include
    -o, --output <output>                    Output file

ARGS:
    <input>    Input VCF file
```

### Replace contig names

```
USAGE:
    sequencetoolkit vcfutils replace-contig [FLAGS] [OPTIONS] --contig-mapping <contig-mapping> [input]

FLAGS:
    -a, --add-chr-prefix       add "chr" prefix to human chromosomes
    -h, --help                 Prints help information
    -r, --remove-chr-prefix    remove "chr" prefix from human chromosomes
    -e, --replace-refseq       Replace RefSeq contig name with chromosome number (for dbSNP 152/153)
    -V, --version              Prints version information

OPTIONS:
    -m, --contig-mapping <contig-mapping>    contig mapping file (csv/tsv)
    -o, --output <output>                    Output file

ARGS:
    <input>    Input VCF file
```

Note: contig mapping file should not include header line. A first column should be
original contig name, and a second column should be new contig name.

### Replace sample names

```
USAGE:
    sequencetoolkit vcfutils replace-sample-names [FLAGS] [OPTIONS] --sample-mapping <mapping> [input]

FLAGS:
    -h, --help                 Prints help information
    -r, --random               Shuffle sample order
    -s, --rename-sequential    Rename sample as "sample_1", "sample_2", ...
    -V, --version              Prints version information

OPTIONS:
    -m, --sample-mapping <mapping>    sample name mapping file (csv/tsv)
    -o, --output <output>             Output VCF file

ARGS:
    <input>    Input VCF file
```

### Rewrite format

Remove unnecessary FORMAT tags.

```
USAGE:
    sequencetoolkit vcfutils rewrite-format [OPTIONS] [--] [input]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -e, --exclude-info-list <exclude>...    Black list of format tags
    -f, --format-list <format>...           White list of format tags
    -o, --output <output>                   Output file

ARGS:
    <input>    Input VCF file
```

### Rewrite INFO

Remove unnecessary INFO tags.

```
USAGE:
    sequencetoolkit vcfutils rewrite-info [OPTIONS] [--] [input]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -e, --exclude-info-list <exclude>...    Black list of format tags
    -i, --info-list <info>...               White list of info tags
    -o, --output <output>                   Output file

ARGS:
    <input>    Input VCF file
```

### Expand BED regions

```
USAGE:
    sequencetoolkit bedutils expand [OPTIONS] --expand <expand> [input]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -e, --expand <expand>                Expand length
    -n, --expand-end <expand-end>        Expand length
    -s, --expand-start <expand-start>    Expand length
    -o, --output <output>                Output Text

ARGS:
    <input>...    Input BED files
```

### Merge overlapped regions

```
USAGE:
    sequencetoolkit bedutils merge [OPTIONS] [input]...

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -o, --output <output>    Output Text

ARGS:
    <input>...    Input BED files
```