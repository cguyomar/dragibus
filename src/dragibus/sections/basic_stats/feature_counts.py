from dragibus.plots import *
from dragibus.stats import *

def subsection_feature_counts(mdFile,genes,transcripts,exons):
    mdFile.new_header(level=2, title='Number of elements of each kind')

    stats_df = basic_stats(genes,transcripts,exons)
    mdFile.write(stats_df.to_markdown(index=True),wrap_width=0)

    mdFile.new_line()

    return(mdFile)