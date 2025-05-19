from dragibus.plots import *
from dragibus.stats import *

def subsection_monoexonic_count(mdFile,transcripts):
    mdFile.new_header(level=2, title='Proportion of monoexonic transcripts')

    monoexonic_df = monoexonic_stats(transcripts)

    mdFile.write(monoexonic_df.to_markdown(index=True, stralign='left',numalign="left"),wrap_width=0)
    mdFile.new_line()
    return(mdFile)