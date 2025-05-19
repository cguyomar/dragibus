# Dragibus : Quality control of annotation files

## Dependencies : 

- bedtools (tested with v2.30.0)
- homer (tested with v4.11)
- pandoc for html output

## Installation

```
git clone https://github.com/cguyomar/dragibus.git
cd dragibus
python setup.py install
# python setup.py install --user # when non root
dragibus -h

```

## Options

`--write_gtf` will write a copy of all input annotation files with a `.dragibus.gtf` extension, with the following modifications:
- New `gene` `transcript` of `intron` features if missing from input file 
- New attributes : 

| Tag                   | Value         | Features   | Description                                                             |
|-----------------------|---------------|------------|-------------------------------------------------------------------------|
| is_canonical          | True of False | intron     |                                                                         |
| all_introns_canonical | True or False | transcript |                                                                         |
| nexons                | Int           | transcript |                                                                         |
| large_internal_exon   | True of False | transcript | Is there a >500nt internal exon? (NA for mono and biexonic transcripts) |
| has_polyA_signal      | True of False | transcript |                                                                         |




## Usage : 

```
dragibus --gtf file1.gtf file2.gtf ... --fasta genome.fasta --out report.md
dragibus --gtf file1.gtf file2.gtf ... --fasta genome.fasta --out report.html
```
