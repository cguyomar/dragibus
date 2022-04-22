import pandas as pd



def basic_stats(genes,transcripts_by_f,exons_by_f):
    
    values = dict()

    for f in exons_by_f.keys():

        exons = exons_by_f[f]
        transcripts = transcripts_by_f[f]

    nb_exons = len(exons)
    nb_distinct_exons = len({(e.chr,e.start,e.end) for e in exons})

    nb_transcripts = len(transcripts)

    nb_genes = len(genes)

    nb_introns = len({i for t in transcripts.values() for i in t.introns})
    nb_distinct_introns = len({(i.chr,i.start,i.end) for t in transcripts.values() for i in t.introns})

        values[f] = [nb_exons,nb_distinct_exons,nb_transcripts,nb_genes,nb_introns,nb_distinct_introns]
    
    
    rownames = ["Number of exons",
                "Number of distinct exons",
                "Number of transcripts",
                "Number of genes",
                "Number of introns",
                "Number of distinct introns"]

    df = pd.DataFrame(values,index=rownames)

    return(df)

def monoexonic_stats(transcripts_by_f):

    files = list(transcripts_by_f.keys())
    values = []
    for f in files:
        transcripts = transcripts_by_f[f]

    nb_monoexonic = len({t for t in transcripts.values() if len(t.introns)==1})
        res = [nb_monoexonic,
                nb_monoexonic/len(transcripts),
                len(transcripts) - nb_monoexonic,
                (len(transcripts) - nb_monoexonic)/len(transcripts)]
        values.append(res)
        # values[f] = values[f] + [i/len(transcripts) for i in values[f]]
        # percentages = [i/len(transcripts) for i in values]

    colnames = ["Monoexonic transcripts","% monoexonic transcripts","Multiexonic transcripts","% multiexonic transcripts"]

    df = pd.DataFrame(values, columns = colnames,index=files)

    # df = pd.DataFrame(list(zip(rownames,values,percentages)),columns = ["Transcript type","Number of transcripts","% of transcripts"])
    df["% monoexonic transcripts"] = df["% monoexonic transcripts"].map("{:.2%}".format)
    df["% multiexonic transcripts"] = df["% multiexonic transcripts"].map("{:.2%}".format)
    return(df)




# def nb_transcripts_by_cdna_length(transcripts,l):
#     filtered_transcripts = {t for t in transcripts.values() if t.cdna_length > l}
#     return(len(filtered_transcripts))



def numeric_feature_distribution(stat_by_fname,breaks):
    # For a set of features ( fname -> features)
    # Compute distribution (nb and perc) of property accessed by fun following breaks

    res = dict()
    res_perc = dict()
    for file in stat_by_fname.keys():
        res[file] = [len({f for f in stat_by_fname[file] if f > l}) for l in breaks]
        res_perc[file] = [r / len(stat_by_fname[file]) for r in res[file]]
    df = pd.DataFrame(res,index=breaks)
    df_perc = pd.DataFrame(res_perc,index=breaks)
    df_perc = df_perc.applymap("{:.2%}".format)

    return(df,df_perc)

def binary_feature_distribution(stat_by_fname,rename={True:"True",False:"False",None:"None"}):
    res = dict()
    for file in stat_by_fname.keys():
        s = pd.value_counts(stat_by_fname[file],dropna=False,sort=True)
        s = s.reindex([True,False,None])
        t = s[True]
        f = s[False]
        row = (t,t/(t+f),f,f/(t+f),s[None])
        res[file] = row

    df = pd.DataFrame.from_dict(res, orient='index',
    columns=(rename[True]+" (nb)",
            rename[True]+" (%)",
            rename[False]+" (nb)",
            rename[False]+" (%)",
            rename[None])
    )
    # Todo : manage key error

    df.iloc[:,1] = df.iloc[:,1].map("{:.2%}".format)
    df.iloc[:,3] = df.iloc[:,3].map("{:.2%}".format)
    return(df)


def collect_stat(features_by_fname,fun,feature_struc="dict"):
    # For a dictionary filename:features, return a dictionary filename:stat computed by fun
    if feature_struc == "dict":
        res = {fname:[fun(f) for f in feature.values()] for fname,feature in features_by_fname.items() }
    elif feature_struc == "set":
        res = {fname:[fun(f) for f in feature] for fname,feature in features_by_fname.items() }
    return(res)




# def nb_transcripts_by_internal_exon_length(transcripts,l):
#     # Number of transcripts with at least one internal exon longer than l
#     # (among transcripts with at least 3 exons)
#     transcripts_3exons = {t for t in transcripts if len(t.exons)>=3}
#     transcripts_3exons_filtered = {t for t in transcripts_3exons if max(e.length for e in t.exons[1:-1]) > l}

#     values = [len(transcripts_3exons_filtered)]
#     perc = [len(transcripts_3exons_filtered)/len(transcripts_3exons)]
#     rownames = ["Transcripts with internal exon(s) longer than "+str(l)]
#     df = pd.DataFrame(list(zip(rownames,values,perc)),columns= ["","Number of transcripts","% of transcripts"])
#     df.iloc[:,2] = df.iloc[:,2].map("{:.2%}".format)

#     return(df)

# def nb_internal_exons_by_length(transcripts,l):
#     # Get number of distinct internal exons longer than l
#     # transcripts_3exons = {t for t in transcripts.values() if len(t.exons)>=3}
#     internal_exons = set()
#     for t in transcripts:
#         internal_exons.update({(e.chr,e.start,e.end) for e in t.exons[1:-1]})
#     nb_internal_exons = len(internal_exons)
#     nb_internal_exons_filtered = len({e for e in internal_exons if (e[2]-e[1])>l})

#     values = [nb_internal_exons_filtered]
#     perc = [nb_internal_exons_filtered / nb_internal_exons]
#     rownames = ["Internal exons longer than "+str(l)]
#     df = pd.DataFrame(list(zip(rownames,values,perc)),columns= ["","Number of unique exons","% of unique exons"])
#     df.iloc[:,2] = df.iloc[:,2].map("{:.2%}".format)

#     return(df)

# def intron_canonicity_stats(transcripts):
#     s = pd.value_counts([i.canonic for t in transcripts for i in t.introns],dropna=False,sort=False)
#     s = s.reindex([True,False,None])
#     df = s.to_frame()
#     df.rename(index={True:"Canonical",False:"Non canonical",None:"NA"},columns={0:"Number of introns"},inplace=True)

#     return(df)
