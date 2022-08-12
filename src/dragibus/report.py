from mdutils.mdutils import MdUtils
from mdutils import Html
import os
import pandas as pd
import plotly
import subprocess
from collections import defaultdict

from dragibus.plots import *
from dragibus.stats import *
from dragibus.sections.basic_stats.basic_stats import *



def make_report(genes,transcripts,exons,introns,errors,mode,skip_polya,out_prefix):

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
    #  Parsing errors
    #
    mdFile.new_header(level=1, title='Parsing errors')

    errors_df = pandas.DataFrame.from_dict(errors)
    mdFile.write(errors_df.to_markdown(index=True),wrap_width=0)
    mdFile.new_line()


    #
    #  Statistics for the whole file
    #
    mdFile = make_basic_stats_section(mdFile,genes,transcripts,exons,mode)

  
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


    mdFile.write(df_intron_canonicity.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)

    if not skip_polya:
    #
    #  PolyA signature
        # 
        mdFile.new_line()
        mdFile.new_header(level=2, title='Polyadenylation sites') 
        mdFile.new_line()

        # All transcripts
        mdFile.write("**Presence of polyA sites for all transcripts** \n ")

        polyA_stat = collect_stat(transcripts, lambda t : t.polyA) 
        df_transcripts_polyA = binary_feature_distribution(polyA_stat,
                rename={True:"PolyA site",
                        False:"No polyA site",
                        None:"NA"}
        )

        mdFile.new_line()
        mdFile.write(df_transcripts_polyA.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
        mdFile.new_line()
        mdFile.new_line()
        
        # Unique tts
        mdFile.write("**Presence of polyA sites for unique transcription terminaison sites** \n ")

        get_tts = lambda t: (t.chr,t.end,t.strand) if t.strand=="+" else (t.chr,t.start,t.strand)
        
        transcripts_unique_tts = defaultdict(dict)
        for f,transcripts_f in transcripts.items():
            unique_tts = { get_tts(t) for t in transcripts_f.values() }
            for id,t in transcripts_f.items():
                tts = get_tts(t)
                if tts in unique_tts:
                    transcripts_unique_tts[f][id] = t
                    unique_tts.remove(tts)
        polyA_stat_unique = collect_stat(transcripts_unique_tts, lambda t : t.polyA) 
        df_transcripts_polyA_unique = binary_feature_distribution(polyA_stat_unique,
                rename={True:"PolyA site",
                        False:"No polyA site",
                        None:"NA"}
        )
        mdFile.new_line()
        mdFile.write(df_transcripts_polyA_unique.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
        mdFile.new_line()
    mdFile.new_line()

    mdFile.create_md_file()

    # make html report from md
    if mode == "html":
        out_name = out_prefix+".html"
        args = [
                'pandoc',
                '-s',
                '-c',
                'pandoc.css',
                '--from',
                'markdown-markdown_in_html_blocks+raw_html',
                '-o',
                out_name,
                out_prefix+".md"
            ]
        subprocess.check_call(args)
