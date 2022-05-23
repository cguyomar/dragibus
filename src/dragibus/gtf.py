from dragibus.features import Feature, Gene, Transcript, Exon



def parse_gtf(in_file):
    # Read a (sorted) gtf file and returns sets of features with a gene -> transcript -> exon structure
    # Transcript/gene structure is inferred from exoon entries
    # Todo : check consistency with gene/transcript lines if present
    
    genes = dict() # gene_id -> gene
    transcripts = dict() # transcript_id -> transcript
    exons = set()

    # Features other than exons preexisting in the input file
    infile_transcripts = dict()
    infile_genes = dict()

    with open(in_file) as input:
        for line in input:
            # print(line)
            if not line.startswith("#"):
                line = line.split("\t")
                assert len(line)==9
                if line[2]=="exon":
                    f = Exon(line[0:8],line[8])
                    exons.add(f)

                    # Add transcript from exon transcript_id                    
                    if f.transcript_id in transcripts.keys():
                        t = transcripts[f.transcript_id]
                        transcripts[f.transcript_id].add_exon(f)
                    else:
                        t = Transcript(line[0:8],{k:v for k,v in f.attributes.items() if k in ["gene_id","transcript_id"]})
                        transcripts[t.id] = t
                        transcripts[f.transcript_id].add_exon(f)

                    # Add gene from exon gene_id                    
                    if f.gene_id not in genes.keys():
                        g = Gene(line[0:8],{"gene_id":f.attributes["gene_id"]})
                        genes[g.id] = g

                    # Associate transcript to gene if not already done
                    if f.transcript_id not in {t.id for t in genes[g.id].transcripts}:
                        genes[t.gene_id].add_transcript(transcripts[f.transcript_id])

                if line[2]=="transcript":
                    f = Transcript(line[0:8],line[8])
                    infile_transcripts[f.id] = f
                    # print(f.attributes)
                if line[2]=="gene":
                    f = Gene(line[0:8],line[8])
                    infile_genes[f.id] = f
                else:
                    f = Feature(line[0:8],line[8])

    # Assign exons to transcripts
    # for e in exons:
    #     assert e.transcript_id in transcripts.keys()
    #     transcripts[e.transcript_id].exons.append(e)

    # Sort exons for each transcript
    for t in transcripts.values():
        t.sort_exons()

    # # Assign transcripts to genes
    # for t in transcripts.values():
    #     assert t.gene_id in genes.keys()
    #     genes[t.gene_id].transcripts.append(t)

    # Compute transcript length
    for t in transcripts.values():
        t.compute_cdna_length()
        t.compute_length()

    # Add intron entries
    for t in transcripts.values():
        t.add_introns()

    return(genes,transcripts,exons)

