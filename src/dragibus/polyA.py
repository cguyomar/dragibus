# import pkgutil
import importlib_resources
import subprocess
import pybedtools
import logging
from pathlib import Path



def scan_genome_for_polyA_motifs(fasta):

    hexamers = []

    pkg = importlib_resources.files("dragibus")
    motif_file = pkg / "assets" / "polyAhex.motif.profile.txt"

    stem = Path(fasta).stem
    outfile_name = stem + '_polya_hexamers.bed'

    if Path(outfile_name).is_file():
        logging.info("Genome has already been searched for hexamers")
    else:
        # with importlib_resources.as_file(motif_file) as motif:
            motif = motif_file
            homer_command = ['scanMotifGenomeWide.pl',motif,fasta,'-keepAll','-bed']
            awk_command = ['awk', "{print $1,$(NF-4)-1,$(NF-3),\".\",\".\",$NF}", "OFS=\t" ]
            f = open(outfile_name, "w")
            p1 = subprocess.Popen(homer_command,stdout=subprocess.PIPE)
            p2 = subprocess.Popen(awk_command,stdin=p1.stdout,stdout=f)
            p2.wait()
            f.close()

    with open(outfile_name) as file:
        for line in file:
            line = line.strip().split("\t")
            hexamers.append([line[0].split(" ")[0],line[-5],line[-4],".",".",line[-1]])

    return(hexamers)

def find_transcripts_with_polya_signal(transcripts,hexamers,flank):
    hexamers_bed = pybedtools.BedTool(hexamers)
    transcripts_to_bed = []
    for t in transcripts.values():
        if t.strand == "-":
            transcripts_to_bed.append((t.chr,max(0,t.start-1),(t.start-1+flank),t.id,".",t.strand))
        elif t.strand == "+":
            transcripts_to_bed.append((t.chr,max(0,t.end-flank),t.end,t.id,".",t.strand))

    transcripts_bed = pybedtools.BedTool(transcripts_to_bed)

    res = transcripts_bed.intersect(hexamers_bed, u=True, s=True)
    transcripts_bed.saveas('a.bed')
    hexamers_bed.saveas('b.bed')
    print(transcripts_bed.count())
    res.saveas("polya.bed")
    print(res.count())
    transcripts_with_polyA = set()
    for rec in res:
        transcripts_with_polyA.add(rec[3])
    for t_id in transcripts.keys():
        transcripts[t_id].polyA = t_id in transcripts_with_polyA

