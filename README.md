# CircosAlignmentPlotter
Converts a part of an alignment (.PAF perhaps others sometimes) to a Circos image using BED and fasta files.

### Requirements
- biopython
- circos
- pandas
- numpy

### Usage
```bash
usage: p2c.py [-h] [--templates TEMPLATES] [--min-length MIN_LENGTH]
              [--bed BED] [--cov COV] [--snps SNPS]
              PAF Query Reference BED Output

Circos plot all alignments against each target in the .paf file.

positional arguments:
  PAF                   An alignment .paf formatted
  Query                 A fasta file containing the query sequences
  Reference             A fasta file containing the target sequences
  BED                   A .bed file containing the target contigs regions to plot
  Output                A directory path

optional arguments:
  -h, --help                              show this help message and exit
  --templates TEMPLATES, -t TEMPLATES     A directory path containing the template files
  --min-length MIN_LENGTH, -m MIN_LENGTH  Minimum length. Default: 1000.
  --bed BED, -b BED                       A .bed file with regions to highlight.
  --scov COV, -sc COV                     A file formatted like the output of samtools depth.
  --ccov COV, -cc COV                     A file formatted for circos.
  --snps SNPS, -s SNPS                    A .vcf file of the query and/or target.
```
