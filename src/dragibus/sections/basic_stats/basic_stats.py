from mdutils.mdutils import MdUtils
from mdutils import Html
import plotly
import pandas as pd


from dragibus.plots import *
from dragibus.stats import *

from dragibus.sections.basic_stats.feature_counts import *
from dragibus.sections.basic_stats.monoexonic_counts import *
from dragibus.sections.basic_stats.transcripts_per_gene import *
from dragibus.sections.basic_stats.exons_per_transcript import *
from dragibus.sections.basic_stats.transcripts_length_distribution import *



def make_basic_stats_section(mdFile,genes,transcripts,exons,mode):

    mdFile.new_header(level=1, title='Quality summary for the whole file')

    mdFile = subsection_feature_counts(mdFile,genes,transcripts,exons)
    mdFile = subsection_monoexonic_count(mdFile,transcripts)
    mdFile = subsection_transcripts_per_gene(mdFile,genes,mode)
    mdFile = subsection_exons_per_transcript(mdFile,transcripts,mode)
    mdFile = subsection_transcripts_length_distribution(mdFile,transcripts,mode)

    return(mdFile)