from mdutils.mdutils import MdUtils
from mdutils import Html
import os
import pandas as pd
import plotly
import subprocess

from dragibus.plots import *
from dragibus.stats import *



def make_report(genes,transcripts,exons,introns,mode,out_prefix):

    if not os.path.exists('ressources'):
        os.mkdir("ressources")


    #
    #  Create markdown
    #
    mdFile = MdUtils(file_name=out_prefix+".md", title='Dragibus evaluation report')


    # Import plotly dependency
    if mode=="html":
        mdFile.write('<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>')

    #
    #  Statistics for the whole file
    #

    mdFile.new_header(level=1, title='Quality summary for the whole file')

    mdFile.new_header(level=2, title='Number of elements of each kind')

    stats_df = basic_stats(genes,transcripts,exons)
    mdFile.write(stats_df.to_markdown(index=True),wrap_width=0)

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
    length_breaks = [50000,100000,500000,1000000,2000000]

    transcript_length_stat = collect_stat(transcripts,lambda t:t.transcript_length)
    
    df_transcript_length, df_transcript_length_perc = numeric_feature_distribution(transcript_length_stat,length_breaks)

    # df_transcript_length = pd.DataFrame.from_dict({l:nb_transcripts_by_transcript_length(transcripts,l) for l in length_breaks},
    #     orient='index',
    #     columns=["Number of transcripts"])
    # df_transcript_length['Percentage of transcripts'] = df_transcript_length["Number of transcripts"] / len(transcripts)
    # df_transcript_length['Percentage of transcripts'] = df_transcript_length['Percentage of transcripts'].map("{:.2%}".format)
    
    mdFile.new_line("\n")
    mdFile.write("**Distribution by number of transcripts** ")
    mdFile.new_line("\n")

    mdFile.write(df_transcript_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    mdFile.new_line()

    mdFile.new_line("\n")
    mdFile.write("**Distribution by percentage of transcripts**")
    mdFile.new_line("\n")

    mdFile.write(df_transcript_length_perc.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()

    #  Non interactive plot - for pdf/md output
    # plot_transcript_length_histogram(transcripts.values(),"ressources/transcript_length_histogram.svg")

    fig_transcript_size = plot_transcript_length_density(transcript_length_stat,5000,"Distribution of transcript size")
    if mode == "html":
        html_transcript_size = plotly.offline.plot(fig_transcript_size, include_plotlyjs=False, output_type='div')
        mdFile.new_paragraph(Html.paragraph(text=html_transcript_size),wrap_width=0)
    elif mode == "markdown":
        fig_transcript_size.write_image("ressources/transcript_length_histogram.svg")
        mdFile.new_paragraph(Html.image(path="ressources/transcript_length_histogram.svg"))
    mdFile.new_line()

    mdFile.new_header(level=3, title='Transcripts cdna length distribution')

    cdna_length_stat = collect_stat(transcripts,lambda t:t.cdna_length)
    df_cdna_length, df_cdna_length_perc = numeric_feature_distribution(cdna_length_stat,[2000,5000,10000,50000,100000])

    mdFile.new_line("\n")
    mdFile.write("**Distribution by number of transcripts** ")
    mdFile.new_line("\n")

    mdFile.write(df_cdna_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    mdFile.new_line()

    mdFile.new_line("\n")
    mdFile.write("**Distribution by percentage of transcripts**")
    mdFile.new_line("\n")

    mdFile.write(df_cdna_length_perc.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    # mdFile.new_line()

    fig_cdna_size = plot_transcript_length_density(cdna_length_stat,500,"Distribution of cdna length size")
    html_cdna_size = plotly.offline.plot(fig_cdna_size, include_plotlyjs=False, output_type='div')
    # fig.write_json("plot.json")
    mdFile.new_paragraph(Html.paragraph(text=html_cdna_size),wrap_width=0)
    mdFile.new_line()

    # Static plot :
    # plot_cdna_length_histogram(transcripts.values(),"ressources/cdna_length_histogram.svg")
    # mdFile.new_paragraph(Html.image(path="ressources/cdna_length_histogram.svg"))
    # mdFile.new_line()


    #
    #  Outlier transcripts
    #
    mdFile.new_header(level=2, title='Outlier transcripts')

    thr_transcript_length = 100000
    thr_internal_exon_length = 500

    mdFile.new_header(level=3, title='Transcripts with at least one internal exon larger than ' + str(thr_internal_exon_length) + " bp.")
    mdFile.new_line()

    large_internal_exon_stat_by_transcript = collect_stat(transcripts,lambda t : None if len(t.exons) < 3  else max(e.length for e in t.exons[1:-1])>thr_internal_exon_length)
    
    df_transcripts_by_internal_exon_length = binary_feature_distribution(large_internal_exon_stat_by_transcript,
            rename={True:"Larger than " + str(thr_internal_exon_length) + " bp.",
                    False:"Smaller than " + str(thr_internal_exon_length) + " bp.",
                    None:"Less than 3 exons"}
    )
   
    # df_transcripts_by_internal_exon_stats = nb_transcripts_by_internal_exon_length(transcripts.values(),thr_internal_exon_length )
    mdFile.write(df_transcripts_by_internal_exon_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line("\n")


    large_internal_exons_stat = collect_stat(exons, lambda e : e.length > thr_internal_exon_length if e.is_internal else None,"set") 
    df_internal_exons_by_length = binary_feature_distribution(large_internal_exons_stat,
            rename={True:"Larger than " + str(thr_internal_exon_length) + " bp.",
                    False:"Smaller than " + str(thr_internal_exon_length) + " bp.",
                    None:"Not internal"}
    )

    mdFile.new_header(level=3, title='Internal exons larger than ' + str(thr_internal_exon_length) + " bp.")
    mdFile.new_line()

    # df_internal_exons_by_length = nb_internal_exons_by_length(transcripts.values(),thr_internal_exon_length)
    mdFile.write(df_internal_exons_by_length.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()


    mdFile.new_header(level=2, title='Intron canonicity')
  
    mdFile.new_line()
    intron_canonicity_stat = collect_stat(introns, lambda i : i.canonic,"set") 
    df_intron_canonicity = binary_feature_distribution(intron_canonicity_stat,
            rename={True:"Canonical ",
                    False:"Non canonical",
                    None:"NA"}
    )


    # df_canonicity = intron_canonicity_stats(transcripts.values())
    mdFile.write(df_intron_canonicity.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)


    mdFile.create_md_file()

    # make html report from md
    if mode == "html":
        out_name = out_prefix+".html"
        args = ['pandoc','-s','-c','pandoc.css','--from','markdown-markdown_in_html_blocks+raw_html','-o','out.html','out.md']
        # subprocess.Popen(args)
        subprocess.check_call(args)
