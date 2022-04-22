#!/usr/bin/env python3

import argparse
import os

import dragibus
from dragibus.stats import collect_stat

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("gtf", help="Input annotation file")
    parser.add_argument("fasta", help="Genome fasta file")
    parser.add_argument("out", help="Output file")
    parser.add_argument("--mode", dest='mode', help="Output format : markdown or html")
    # parser.add_argument("--out",dest='out_file')
    # parser.add_argument("--fasta",dest='fasta')

    args = parser.parse_args()

    # in_file = "small.gtf"
    # fasta = "sus_scrofa.fa"
    # in_file = "novel.gtf"
    in_file = args.gtf
    fasta = args.fasta    
    out_file = args.out
    out_prefix,ext = os.path.splitext(out_file)

    if args.mode:
        mode = args.mode

    if ext==".md":
        mode = "markdown"
    elif ext==".html":
        mode = "html"
        
    if mode not in ["html","markdown"]:
        print("Please choose a valid value for --mode : html or markdown")
        exit(1)
    if mode=="markdown":
        ext=".md"
    elif mode=="html":
        ext=".html"

    

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


    dragibus.make_report(genes,transcripts,exons,introns,mode,out_prefix)



if __name__ == "__main__":
    main()





        
