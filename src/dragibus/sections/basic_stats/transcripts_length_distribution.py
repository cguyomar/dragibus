import plotly
from mdutils.mdutils import MdUtils
from mdutils import Html

from dragibus.plots import *
from dragibus.stats import *

def subsection_transcripts_length_distribution(mdFile,transcripts,mode):
   
    mdFile.new_header(level=2, title='Transcripts length distribution')

    #
    # transcript length distribution
    #
    mdFile.new_header(level=3, title='Transcripts size distribution')
    length_breaks = [50000,100000,500000,1000000,2000000]

    transcript_length_stat = collect_stat(transcripts,lambda t:t.transcript_length)
    df_transcript_length = numeric_feature_distribution(transcript_length_stat,length_breaks)

    mdFile.write("Number of transcripts longer than n nt.")
    mdFile.new_line("\n")
    mdFile.write(df_transcript_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    mdFile.new_line()

    #  Non interactive plot - for pdf/md output
    # plot_transcript_length_histogram(transcripts.values(),"ressources/transcript_length_histogram.svg")
    up_limit = max([np.quantile(transcript_length_stat[f],0.95) for f in transcript_length_stat.keys()])
    # bin_width = int(up_limit/30)
    # fig_transcript_size = plot_transcript_length_density(transcript_length_stat,bin_width,up_limit,"Distribution of transcript size")
    df = feature_long_table(transcript_length_stat)
    fig_transcript_size = plot_histogram_distribution(df,max=up_limit,title="Distribution of transcript size",log=False)

    if mode == "html":
        html_transcript_size = plotly.offline.plot(fig_transcript_size, include_plotlyjs=False, output_type='div')
        mdFile.new_paragraph(Html.paragraph(text=html_transcript_size),wrap_width=0)
    elif mode == "markdown":
        fig_transcript_size.write_image("ressources/transcript_length_histogram.svg")
        mdFile.new_paragraph(Html.image(path="ressources/transcript_length_histogram.svg"))
    mdFile.new_line()

    #
    # cDNA length distribution
    #
    mdFile.new_header(level=3, title='Transcripts cdna length distribution')

    cdna_length_stat = collect_stat(transcripts,lambda t:t.cdna_length)
    up_limit_cdna_length = max([np.quantile(cdna_length_stat[f],0.95) for f in cdna_length_stat.keys()])
    # bin_width = int(up_limit/30)
    
    df_cdna_length = numeric_feature_distribution(cdna_length_stat,[2000,5000,10000,50000,100000])

    mdFile.new_line("\n")
    mdFile.write("**Distribution of cdna length** ")
    mdFile.new_line("\n")

    mdFile.write("Number of cdnas longer than n nt. \n")
    mdFile.new_line("\n")
    mdFile.write(df_cdna_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    mdFile.new_line()

    # fig_cdna_size = plot_transcript_length_density(cdna_length_stat,bin_width,up_limit,"Distribution of cdna length size")
    df = feature_long_table(cdna_length_stat)
    fig_cdna_size = plot_histogram_distribution(df,max=up_limit_cdna_length,title="Distribution of cdna length size",log=False)

    if mode == "html":
        html_cdna_size = plotly.offline.plot(fig_cdna_size, include_plotlyjs=False, output_type='div')
        mdFile.new_paragraph(Html.paragraph(text=html_cdna_size),wrap_width=0)
    elif mode == "markdown":
        fig_cdna_size.write_image


    return(mdFile)