import pandas as pd
from collections import Counter
import numpy as np
from itertools import chain



def basic_stats(genes_by_f,transcripts_by_f,exons_by_f):
    
    values = dict()

    for f in exons_by_f.keys():

        exons = exons_by_f[f]
        transcripts = transcripts_by_f[f]
        genes = genes_by_f[f]

        nb_CDS = len([t for t in transcripts.values() if t.CDS != None])

        nb_exons = len(exons)
        nb_distinct_exons = len({(e.chr,e.start,e.end) for e in exons})
        nb_transcripts = len(transcripts)
        nb_genes = len(genes)
        nb_introns = len({i for t in transcripts.values() for i in t.introns})
        nb_distinct_introns = len({(i.chr,i.start,i.end) for t in transcripts.values() for i in t.introns})

        values[f] = [
            nb_genes,
            nb_transcripts,
            nb_CDS,
            nb_exons,
            nb_distinct_exons,
            nb_introns,
            nb_distinct_introns
            ]
    
    
    rownames = [
        "Number of genes",
        "Number of transcripts",
        "Number of CDS",
        "Number of exons",
        "Number of distinct exons",
        "Number of introns",
        "Number of distinct introns"
        ]

    df = pd.DataFrame(values,index=rownames)

    return(df)

def monoexonic_stats(transcripts_by_f):

    files = list(transcripts_by_f.keys())
    values = []
    for f in files:
        transcripts = transcripts_by_f[f]

        nb_monoexonic = len({t for t in transcripts.values() if len(t.exons)==1})
        res = [nb_monoexonic,
                nb_monoexonic/len(transcripts),
                len(transcripts) - nb_monoexonic,
                (len(transcripts) - nb_monoexonic)/len(transcripts)]
        values.append(res)

    colnames = ["Monoexonic transcripts","% monoexonic transcripts","Multiexonic transcripts","% multiexonic transcripts"]

    df = pd.DataFrame(values, columns = colnames,index=files)

    df["% monoexonic transcripts"] = df["% monoexonic transcripts"].map("{:.2%}".format)
    df["% multiexonic transcripts"] = df["% multiexonic transcripts"].map("{:.2%}".format)
    return(df)

def discrete_feature_distribution(stat_by_fname):

    df = pd.DataFrame()
    for file in stat_by_fname.keys():
        count = Counter(stat_by_fname[file])
        temp = pd.DataFrame.from_dict(count, orient='index')
        df = pd.merge(df,temp,left_index=True,right_index=True,how="outer")
        # Name column
        df.set_axis([*df.columns[:-1], file], axis=1, inplace=True)
    df.sort_index(inplace=True)
    return(df)

def feature_long_table(stat_by_fname):
    # Return long format tab
    df = pd.DataFrame({
        'file' : np.repeat(list(stat_by_fname.keys()), [len(v) for k, v in stat_by_fname.items()]),
        'value' : list(chain.from_iterable(stat_by_fname.values())), 
    })
    return(df)


def numeric_feature_distribution(stat_by_fname,breaks):
    # For a set of features ( fname -> features)
    # Compute distribution (nb and perc) of property accessed by fun following breaks

    res = dict()
    # res_perc = dict()
    for file in stat_by_fname.keys():
        res[file] = [len([f for f in stat_by_fname[file] if f >= l]) for l in breaks]
        res[file+" (%)"] = map("{:.2%}".format,[r / len(stat_by_fname[file]) for r in res[file]])
    df = pd.DataFrame(res,index=breaks)
    # df_perc = pd.DataFrame(res_perc,index=breaks)
    # df_perc = df_perc.applymap("{:.2%}".format)
    # print(df)

    return(df)

def binary_feature_distribution(stat_by_fname,rename={True:"True",False:"False",None:"None"}):
    res = dict()
    for file in stat_by_fname.keys():
        s = pd.Series(stat_by_fname[file]).value_counts(dropna=False, sort=True)
        s = s.reindex([True,False,None], fill_value=0)
        t = s.get(True, 0)
        f = s.get(False, 0)
        row = (
            t,
            "{:.2%}".format(t/(t+f)),
            f,
            "{:.2%}".format(f/(t+f)),
            s[None])
        res[file] = row

    df = pd.DataFrame.from_dict(res, orient='index',
    columns=(rename[True]+" (nb)",
            rename[True]+" (%)",
            rename[False]+" (nb)",
            rename[False]+" (%)",
            rename[None])
    )
    # Drop NA columns if not found
    if sum(df[rename[None]].notna())==0:
        df.drop(axis=1,labels=rename[None],inplace=True)

    return(df)


def collect_stat(features_by_fname,fun,feature_struc="dict"):
    # For a dictionary filename:features, return a dictionary filename:stat computed by fun
    if feature_struc == "dict":
        res = {fname:[fun(f) for f in feature.values()] for fname,feature in features_by_fname.items() }
    elif feature_struc == "set":
        res = {fname:[fun(f) for f in feature] for fname,feature in features_by_fname.items() }
    return(res)


