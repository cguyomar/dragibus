from dragibus.plots import *
from dragibus.stats import *

def subsection_CDS_counts(mdFile,genes,transcripts,exons):
    mdFile.new_header(level=2, title='Number of transcripts with CDS')

    hasCDS_by_transcript = collect_stat(transcripts,lambda t : t.CDS )
    
    df_transcripts_with_CDS = binary_feature_distribution(hasCDS_by_transcript,
            rename={True:"Transcripts with CDS",
                    False:"Transcripts without CDS",
                    None:"NA"}
    )

    mdFile.write(df_transcripts_with_CDS.to_markdown(index=True),wrap_width=0)

    mdFile.new_line()

    return(mdFile)