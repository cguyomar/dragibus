# import pkgutil
import importlib_resources
import subprocess
import pybedtools
from pathlib import Path



def scan_genome_for_polyA_motifs(fasta):

    hexamers = []

    pkg = importlib_resources.files("dragibus")
    motif_file = pkg / "assets" / "polyAhex.motif.profile.txt"

    stem = Path(fasta).stem
    outfile_name = stem + '_polya_hexamers.bed'

    if Path(outfile_name).is_file():
        print("Genome has already been searched for hexamers")
    else:
        homer_command = ['scanMotifGenomeWide.pl',motif_file,fasta,'-keepAll','-bed']
        awk_command = ['awk', "{print $1,$(NF-4)-1,$(NF-3),\".\",\".\",$NF}", "OFS='\t'" ]
        f = open(outfile_name, "w")
        subprocess.run(homer_command,stdout=f)
        print(awk_command)

    with open(outfile_name) as file:
        for line in file:
            # line = line.split("\t")
            line = line.strip().split("\t")
            # print(line)

            hexamers.append([line[0].split(" ")[0],line[-5],line[-4],".",".",line[-1]])

            # f.write("\t".join([line[0],line[-5],line[-4],".",".",line[-1]])+"\n")
            # print $1,$(NF-4)-1,$(NF-3),".",".",$NF

    pybedtools.BedTool(hexamers).saveas("toto.bed")

    return(hexamers)

def find_transcripts_with_polya_signal(transcripts,hexamers,flank):
    hexamers_bed = pybedtools.BedTool(hexamers)
    transcripts_to_bed = []
    for t in transcripts.values():
        if t.strand == "-":
            transcripts_to_bed.append((t.chr,t.start-1-flank,t.end,t.id,".",t.strand))
            print(t.start-1-flank,t.end,t.id,".",t.strand)
        elif t.strand == "+":
            transcripts_to_bed.append((t.chr,t.start-1,t.end+flank,t.id,".",t.strand))
            print(t.start-1,t.end+flank,t.id,".",t.strand)
    transcripts_bed = pybedtools.BedTool(transcripts_to_bed)

    res = transcripts_bed.intersect(hexamers_bed, u=True, s=True)
    transcripts_with_polyA = set()
    for rec in res:
        transcripts_with_polyA.add(rec[3])
    for t_id in transcripts.keys():
        transcripts[t_id].polyA = t_id in transcripts_with_polyA

