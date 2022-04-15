import pybedtools

def get_sequence(chr,start,end,fasta):
    if start < end:
        bed_entry = "\t".join([chr,str(start-1),str(end)])
        bed = pybedtools.BedTool(bed_entry,from_string=True)
        fa = bed.sequence(fi=fasta,s=True).seqfn
        seq = open(fa).read().split("\n")[1]
        return(seq)
    else:
        bed_entry = "\t".join([chr,str(end-1),str(start)])
        bed = pybedtools.BedTool(bed_entry,from_string=True)
        fa = bed.sequence(fi=fasta).seqfn
        seq = open(fa).read().split("\n")[1]
        seq = str(Seq(seq).reverse_complement())
        return(seq)