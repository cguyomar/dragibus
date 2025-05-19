import plotly
from mdutils.mdutils import MdUtils
from mdutils import Html

from dragibus.plots import *
from dragibus.stats import *

def subsection_exons_per_transcript(mdFile,transcripts,mode):
   
    mdFile.new_header(level=2, title='Number of exons per transcript')
    mdFile.new_line()

    nb_e_per_transcript = collect_stat(transcripts,lambda t:len(t.exons))
    df_nb_e_per_t = feature_long_table(nb_e_per_transcript)
    up_limit = max([np.quantile(nb_e_per_transcript[f],0.95) for f in nb_e_per_transcript.keys()])
    plot_nb_e_per_t = plot_violin_distribution(df_nb_e_per_t,up_limit)
    if mode == "html":
        html_nb_e_per_t = plotly.offline.plot(plot_nb_e_per_t, include_plotlyjs=False, output_type='div')
        mdFile.new_paragraph(Html.paragraph(text=html_nb_e_per_t),wrap_width=0)
    elif mode == "markdown":
        plot_nb_e_per_t.write_image("ressources/nb_exons_per_transcript.svg")
        mdFile.new_paragraph(Html.image(path="ressources/nb_exons_per_transcript.svg"))
    mdFile.new_line()
    return(mdFile)