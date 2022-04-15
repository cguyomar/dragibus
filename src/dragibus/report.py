from mdutils.mdutils import MdUtils
from mdutils import Html
import os
import pandas as pd

from dragibus.plots import *
from dragibus.stats import *


def make_report(genes,transcripts,exons,filename="out.md"):

    if not os.path.exists('ressources'):
        os.mkdir("ressources")


    #
    #  Create markdown
    #
    mdFile = MdUtils(file_name=filename, title='GFF evaluation report')

    mdFile.new_header(level=1, title='Parsing report')

    #
    #  Statistics for the whole file
    #

    mdFile.new_header(level=1, title='Quality summary for the whole file')

    mdFile.new_header(level=2, title='Number of elements of each kind')

    stats_df = basic_stats(genes,transcripts,exons)
    mdFile.write(stats_df.to_markdown(index=False),wrap_width=0)

    mdFile.new_line()
    mdFile.new_line()

    mdFile.new_header(level=2, title='Proportion of monoexonic transcripts')

    monoexonic_df = monoexonic_stats(transcripts)

    mdFile.write(monoexonic_df.to_markdown(index=False, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()

    #
    #  Length distribution
    #

    mdFile.new_header(level=2, title='Transcripts length distribution')


    mdFile.new_header(level=3, title='Transcripts size distribution')

    df_transcript_length = pd.DataFrame.from_dict({l:nb_transcripts_by_transcript_length(transcripts,l) for l in [50000,100000,500000,1000000,2000000]},
        orient='index',
        columns=["Number of transcripts"])
    df_transcript_length['Percentage of transcripts'] = df_transcript_length["Number of transcripts"] / len(transcripts)
    df_transcript_length['Percentage of transcripts'] = df_transcript_length['Percentage of transcripts'].map("{:.2%}".format)

    mdFile.write(df_transcript_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()

    plot_transcript_length_histogram(transcripts.values(),"ressources/transcript_length_histogram.svg")
    mdFile.new_paragraph(Html.image(path="ressources/transcript_length_histogram.svg"))
    mdFile.new_line()

    mdFile.new_header(level=3, title='Transcripts cdna length distribution')

    df_cdna_length = pd.DataFrame.from_dict({l:nb_transcripts_by_cdna_length(transcripts,l) for l in [2000,5000,10000,50000,100000]},
        orient='index',
        columns=["Number of transcripts"])
    df_cdna_length['Percentage of transcripts'] = df_cdna_length["Number of transcripts"] / len(transcripts)
    df_cdna_length['Percentage of transcripts'] = df_cdna_length['Percentage of transcripts'].map("{:.2%}".format)


    mdFile.write(df_cdna_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()

    plot_cdna_length_histogram(transcripts.values(),"ressources/cdna_length_histogram.svg")
    mdFile.new_paragraph(Html.image(path="ressources/cdna_length_histogram.svg"))
    mdFile.new_line()


    #
    #  Outlier transcripts
    #
    mdFile.new_header(level=2, title='Outlier transcripts')

    thr_transcript_length = 100000
    thr_internal_exon_length = 500

    df_transcripts_by_internal_exon_stats = nb_transcripts_by_internal_exon_length(transcripts.values(),thr_internal_exon_length )
    mdFile.write(df_transcripts_by_internal_exon_stats.to_markdown(index=False, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    mdFile.new_line()


    df_internal_exons_by_length = nb_internal_exons_by_length(transcripts.values(),thr_internal_exon_length)
    mdFile.write(df_internal_exons_by_length.to_markdown(index=False, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()


    mdFile.new_header(level=2, title='Intron canonicity')
    df_canonicity = intron_canonicity_stats(transcripts.values())
    mdFile.write(df_canonicity.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)


    mdFile.create_md_file()