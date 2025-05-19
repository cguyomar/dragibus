import plotly
from mdutils.mdutils import MdUtils
from mdutils import Html

from dragibus.plots import *
from dragibus.stats import *

def subsection_transcripts_per_gene(mdFile,genes,mode):
   
    mdFile.new_header(level=2, title='Number of transcripts per gene')
    mdFile.new_line()

    nb_t_per_gene = collect_stat(genes,lambda g:len(g.transcripts))
    df_nb_t_per_gene = feature_long_table(nb_t_per_gene)
    
    # Find 95% quantile to limit plot
    up_limit = max([np.quantile(nb_t_per_gene[f],0.95) for f in nb_t_per_gene.keys()])

    plot_nb_t_per_gene =  plot_violin_distribution(df_nb_t_per_gene,max=up_limit)
    if mode == "html":
        html_nb_t_per_gene = plotly.offline.plot(plot_nb_t_per_gene, include_plotlyjs=False, output_type='div')
        mdFile.new_paragraph(Html.paragraph(text=html_nb_t_per_gene),wrap_width=0)
    elif mode == "markdown":
        plot_nb_t_per_gene.write_image("ressources/nb_transcripts_per_gene.svg")
        mdFile.new_paragraph(Html.image(path="ressources/nb_transcripts_per_gene.svg"))
    mdFile.new_line()


    return(mdFile)