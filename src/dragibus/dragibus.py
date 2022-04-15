#!/usr/bin/env python3

import argparse

import dragibus

def main():

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--in",dest='in_file')
    parser.add_argument("--out",dest='out_file')
    parser.add_argument("--fasta",dest='fasta')

    # in_file = "small.gtf"
    # fasta = "sus_scrofa.fa"
    # in_file = "novel.gtf"

    genes,transcripts,exons = dragibus.parse_gtf(in_file)

    dragibus.find_canonic_introns(transcripts,fasta)

    dragibus.make_report(genes,transcripts,exons,out_file)


if __name__ == "__main__":
    main()





        
