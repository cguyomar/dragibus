import pandas as pd



def basic_stats(genes,transcripts,exons):
    nb_exons = len(exons)
    nb_distinct_exons = len({(e.chr,e.start,e.end) for e in exons})

    nb_transcripts = len(transcripts)

    nb_genes = len(genes)

    nb_introns = len({i for t in transcripts.values() for i in t.introns})
    nb_distinct_introns = len({(i.chr,i.start,i.end) for t in transcripts.values() for i in t.introns})

    values = [nb_exons,nb_distinct_exons,nb_transcripts,nb_genes,nb_introns,nb_distinct_introns]
    rownames = ["Number of exons",
                "Number of distinct exons",
                "Number of transcripts",
                "Number of genes",
                "Number of introns",
                "Number of distinct introns"]

    df = pd.DataFrame(list(zip(rownames,values)),columns=['Metric','Number of features'])

    # print("\t".join(["nbex","nbdistinctex","nbtr","nbgn","nbintrons","nbdistinctintrons"]))
    # print("\t".join([str(nb_exons),
    #                 str(nb_distinct_exons),
    #                 str(nb_transcripts),
    #                 str(nb_genes),
    #                 str(nb_introns),
    #                 str(nb_distinct_introns)
    #                 ]))

    return(df)

def monoexonic_stats(transcripts):
    nb_monoexonic = len({t for t in transcripts.values() if len(t.introns)==1})
    values = [nb_monoexonic,len(transcripts) - nb_monoexonic]
    percentages = [i/len(transcripts) for i in values]
    rownames = ["Monoexonic transcripts","Multiexonic transcripts"]

    df = pd.DataFrame(list(zip(rownames,values,percentages)),columns = ["Transcript type","Number of transcripts","% of transcripts"])
    df["% of transcripts"] = df["% of transcripts"].map("{:.2%}".format)

    return(df)




def nb_transcripts_by_cdna_length(transcripts,l):
    filtered_transcripts = {t for t in transcripts.values() if t.cdna_length > l}
    return(len(filtered_transcripts))

def nb_transcripts_by_transcript_length(transcripts,l):
    filtered_transcripts = {t for t in transcripts.values() if t.transcript_length > l}
    return(len(filtered_transcripts))

def nb_transcripts_by_internal_exon_length(transcripts,l):
    # Number of transcripts with at least one internal exon longer than l
    # (among transcripts with at least 3 exons)
    transcripts_3exons = {t for t in transcripts if len(t.exons)>=3}
    transcripts_3exons_filtered = {t for t in transcripts_3exons if max(e.length for e in t.exons[1:-1]) > l}

    values = [len(transcripts_3exons_filtered)]
    perc = [len(transcripts_3exons_filtered)/len(transcripts_3exons)]
    rownames = ["Transcripts with internal exon(s) longer than "+str(l)]
    df = pd.DataFrame(list(zip(rownames,values,perc)),columns= ["","Number of transcripts","% of transcripts"])
    df.iloc[:,2] = df.iloc[:,2].map("{:.2%}".format)

    return(df)

def nb_internal_exons_by_length(transcripts,l):
    # Get number of distinct internal exons longer than l
    # transcripts_3exons = {t for t in transcripts.values() if len(t.exons)>=3}
    internal_exons = set()
    for t in transcripts:
        internal_exons.update({(e.chr,e.start,e.end) for e in t.exons[1:-1]})
    nb_internal_exons = len(internal_exons)
    nb_internal_exons_filtered = len({e for e in internal_exons if (e[2]-e[1])>l})

    values = [nb_internal_exons_filtered]
    perc = [nb_internal_exons_filtered / nb_internal_exons]
    rownames = ["Internal exons longer than "+str(l)]
    df = pd.DataFrame(list(zip(rownames,values,perc)),columns= ["","Number of unique exons","% of unique exons"])
    df.iloc[:,2] = df.iloc[:,2].map("{:.2%}".format)

    return(df)

def intron_canonicity_stats(transcripts):
    s = pd.value_counts([i.canonic for t in transcripts for i in t.introns],dropna=False,sort=False)
    s = s.reindex([True,False,None])
    df = s.to_frame()
    df.rename(index={True:"Canonical",False:"Non canonical",None:"NA"},columns={0:"Number of introns"},inplace=True)

    return(df)
