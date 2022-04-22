#!/usr/bin/env python3

import argparse

import dragibus
from dragibus.stats import collect_stat

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--in",dest='in_file')
    parser.add_argument("--out",dest='out_file')
    parser.add_argument("--fasta",dest='fasta')

    args = parser.parse_args()

    # in_file = "small.gtf"
    # fasta = "sus_scrofa.fa"
    # in_file = "novel.gtf"
    in_file = args.in_file
    fasta = args.fasta    
    out_file = args.out_file

    annotation_files = [in_file]

    genes = dict()
    transcripts = dict()
    exons = dict()
    introns = dict()
    for f in annotation_files:
        introns[f] = set()
        genes[f],transcripts[f],exons[f] = dragibus.parse_gtf(in_file)
        # Enrich transcripts with intron information
        dragibus.find_canonic_introns(transcripts[f],fasta)
        for t in transcripts[f].values():
            for i in t.introns:
                introns[f].add(i)


    dragibus.make_report(genes,transcripts,exons,introns,out_file)



if __name__ == "__main__":
    main()





        
