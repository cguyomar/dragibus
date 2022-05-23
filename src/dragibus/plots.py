import pandas
# import seaborn as sns
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import plotly.express as px

def plot_nb_transcripts_per_gene(genes,outname=""):
    df_g = pandas.DataFrame(
        {"nb_transcripts_per_gene" : [len(g.transcripts) for g in genes.values()]}
    )
    nb_transcripts = seaborn.countplot(x='nb_transcripts_per_gene', data=df_g)
    if outname=="":
        plt.show()
    else:
        plt.savefig(outname)
        plt.cla()

def plot_nb_exons_per_transcript(transcripts,outname=""):
    df_t = pandas.DataFrame(
        {"nb_exons_per_transcript" : [len(t.exons) for t in transcripts.values()]}
    )
    nb_exons = seaborn.countplot(x='nb_exons_per_transcript', data=df_t)
    if outname=="":
        plt.show()
    else:
        plt.savefig(outname)
        plt.cla()

# def plot_transcript_length_histogram(transcripts,outname=""):
#     h = sns.histplot({t.transcript_length for t in transcripts},log_scale=(False,True))
#     h.set_title("Histogram of transcript length distribution")
#     h.set_xlabel("Transcript length")
#     if outname=="":
#         plt.show()
#     else:
#         plt.savefig(outname)
#         plt.cla()

def plot_transcript_length_density(features_stat,bin_size,title=""):
    # input : result of collect_stat function with transcript_length
    fnames, lengths = zip(*features_stat.items())

    fig = ff.create_distplot(lengths, fnames, show_hist=True, show_rug=False, bin_size=bin_size)
    fig.layout.update(title=title)
    return(fig)


def plot_discrete_distribution(df,max=10,title=""):
    fig = px.bar(df,barmode='group')
    fig.layout.update(title=title)
    return(fig)

def plot_cdna_length_histogram(transcripts,outname=""):
    h = sns.histplot({t.cdna_length for t in transcripts})
    h.set_title("Histogram of cdna length distribution")
    h.set_xlabel("Cdna length")
    # h.set_xscale('symlog')
    if outname=="":
        plt.show()
    else:
        plt.savefig(outname)
        plt.cla()
