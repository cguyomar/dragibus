import logging
from collections import defaultdict


from dragibus.features import *


def parse_gtf(in_file,errors):
    # Read a (sorted) gtf file and returns sets of features with a gene -> transcript -> exon structure
    # Transcript/gene structure is inferred from exon entries
    # Todo : check consistency with gene/transcript lines if present
    genes = dict() # gene_id -> gene
    transcripts = dict() # transcript_id -> transcript
    exons = set()

    logging.info('Parsing the gtf')

    errors = defaultdict(int)

    # Features other than exons preexisting in the input file
    infile_transcripts = dict()
    infile_genes = dict()

    with open(in_file) as input:
        for line in input:
            if not line.startswith("#"):
                line = line.split("\t")
                try:
                    if len(line)!=9:
                        errors["Malformed gtf line"] += 1
                        raise ValueError("Incorrect number of fields in file - skipping line")
                except ValueError as ve:
                    print(ve)
                    print(line)

                if line[2]=="exon":
                    try:
                        f = Exon(line[0:8],line[8])
                    except InvalidCoordinatesError as e:
                        logging.warning(e)
                        errors[e.key] += 1
                        continue

                    # Add transcript from exon transcript_id                    
                    if f.transcript_id in transcripts.keys():
                        t = transcripts[f.transcript_id]
                        try:
                            transcripts[f.transcript_id].add_exon(f)
                        except DragibusException as e:
                            logging.warning(e)
                            errors[e.key] += 1
                            
                            continue

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
                        try:
                            genes[t.gene_id].add_transcript(transcripts[f.transcript_id])
                        except WrongChromosomeTranscriptError as e:
                            logging.warning(e)
                            errors[e.key] += 1
                        except WrongStrandTranscriptError as e:
                            logging.warning(e)
                            errors[e.key] += 1

                    exons.add(f)

                if line[2]=="transcript":
                    try:
                        f = Transcript(line[0:8],line[8])
                    except DragibusException as e:
                        logging.warning(e)
                        errors[e.key] += 1
                        continue
                    infile_transcripts[f.id] = f
                    # print(f.attributes)
                if line[2]=="gene":
                    try:
                        f = Gene(line[0:8],line[8])
                    except DragibusException as e:
                        logging.warning(e)
                        errors[e.key] += 1
                        continue
                    infile_genes[f.id] = f
                else:
                    try:
                        f = Feature(line[0:8],line[8])
                    except InvalidCoordinatesError as e:
                        print("Invalid")

    # Assign exons to transcripts
    # for e in exons:
    #     assert e.transcript_id in transcripts.keys()
    #     transcripts[e.transcript_id].exons.append(e)

    # Sort exons for each transcript
    for t in transcripts.values():
        try:
            t.sort_exons()
        except DragibusException as e:
            logging.warning(e)
            errors[e.key] += 1
            continue

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
        try:
            t.add_introns()
        except Exception as e:
            print(e)

    # Find and remove transcripts/exons that are orphans due to syntax checking
    transcripts = {t_id:t for t_id,t in transcripts.items() if t.gene_id in genes.keys()}
    exons = {e for e in exons if e.gene_id in genes.keys() and e.transcript_id in transcripts.keys()}

    return(genes,transcripts,exons,errors)

