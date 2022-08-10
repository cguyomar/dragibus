#!/usr/bin/env python3

import argparse
import os
import logging
from collections import defaultdict

import dragibus
from dragibus.stats import collect_stat

def main():

    logging.basicConfig(level=logging.INFO)
    logging.info('Starting Dragibus')

    parser = argparse.ArgumentParser(description='Quality control of annotation files')
    parser.add_argument('--gtf', nargs='+', help='Input annotation file', required=True)
    parser.add_argument("--fasta", help="Genome fasta file",required=True)
    parser.add_argument("--out", help="Output file",required=True)
    parser.add_argument("--mode", dest='mode', help="Output format : markdown or html")
    # parser.add_argument("--out",dest='out_file')
    # parser.add_argument("--fasta",dest='fasta')

    skip_polya=False

    args = parser.parse_args()

    # in_file = "small.gtf"
    # fasta = "sus_scrofa.fa"
    # in_file = "novel.gtf"
    annotation_files = args.gtf
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


    if not skip_polya:
        hexamers = dragibus.scan_genome_for_polyA_motifs(fasta)
   
    errors = defaultdict(int)
    genes = dict()
    transcripts = dict()
    exons = dict()
    introns = dict()
    for f in annotation_files:
        introns[f] = set()
        genes[f],transcripts[f],exons[f],errors = dragibus.parse_gtf(f,errors)
        # Enrich transcripts with intron information
        dragibus.find_canonic_introns(transcripts[f],fasta)
        for t in transcripts[f].values():
            for i in t.introns:
                introns[f].add(i)

    if not skip_polya:
        for f in annotation_files:
            dragibus.find_transcripts_with_polya_signal(transcripts[f],hexamers,10)
    
    dragibus.make_report(genes,transcripts,exons,introns,errors,mode,skip_polya,out_prefix)



if __name__ == "__main__":
    main()





        
