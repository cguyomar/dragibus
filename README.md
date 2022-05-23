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

## Usage : 

```
dragibus --gtf file1.gtf file2.gtf ... --fasta genome.fasta --out report.md
dragibus --gtf file1.gtf file2.gtf ... --fasta genome.fasta --out report.html
```
