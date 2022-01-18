## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo=TRUE, warning=FALSE, message=FALSE, fig.pos = 'h'
)

## ---- fig.show='hold'---------------------------------------------------------
library(diffSeqPatterns)

## ----install_proBatch, fig.show='hold', eval = FALSE--------------------------
#  #install the development version from GitHub:
#  install.packages('devtools')
#  devtools::install_github('ChloeHJ/diffSeqPatterns', build_vignettes = TRUE)

## ---- fig.show='hold'---------------------------------------------------------
library(diffSeqPatterns)
data('DMF5_pos_peptides', 'DMF5_neg_peptides', 'DMF5_antigen_table',  package = 'diffSeqPatterns')

## ---- fig.show='hold'---------------------------------------------------------
diff_pssm_mtx <- compute_diff_pssm(DMF5_pos_peptides, DMF5_neg_peptides)
knitr::kable(diff_pssm_mtx)

## ---- fig.show='hold', fig.width=4, fig.height=2.5----------------------------
plot_raster_diff_pssm(diff_pssm_mtx)
plot_seqlogo_diff_pssm(diff_pssm_mtx)

## ---- fig.show='hold', fig.width=10, fig.height=3-----------------------------
ngram_df <- compute_ngrams(DMF5_pos_peptides, DMF5_neg_peptides, ngram_lengths = c(2, 3, 4, 5))
knitr::kable(head(ngram_df, 5))

## ---- fig.show='hold', fig.width=4, fig.height=2.5----------------------------
plot_point_ngrams(ngram_df, ratio_threshold = 10, n_threshold = 10) 

## ---- fig.show='hold', fig.width=5, fig.height=2------------------------------
plot_bar_top_ngrams(ngram_df, top_n = 20)

## ---- fig.show='hold'---------------------------------------------------------
positional_kmer_data <- compute_positional_kmers(DMF5_pos_peptides, DMF5_neg_peptides)
knitr::kable(head(positional_kmer_data, 5))


## ---- fig.show='hold', fig.width=4, fig.height=2.5----------------------------
plot_point_positional_kmers(positional_kmer_data, ratio_threshold = 25,  n_threshold = 10)

## ---- fig.show='hold', fig.width=7.5, fig.height=2----------------------------
plot_bar_posiitonal_kmers(positional_kmer_data, top_n = 20) 

## ---- fig.show='hold'---------------------------------------------------------
alignment_matrix <- compute_pairwise_alignment(peptides = c(DMF5_pos_peptides, DMF5_neg_peptides[1:20]))
distance_mtx <- compute_pairwise_distance(peptides = c(DMF5_pos_peptides))

knitr::kable(alignment_matrix[1:8, 1:8])


## ---- fig.show='hold', fig.width=8, fig.height=7------------------------------
library(tibble)
library(dplyr)
col_data <- DMF5_antigen_table %>% column_to_rownames(var = 'ContactPositions' )
plot_heatmap_alignment_mtx(alignment_matrix, col_data)
#plot_heatmap_distance_mtx(distance_mtx, col_data)

## ---- fig.show='hold', fig.width=6, fig.height=5------------------------------
plot_network_alignment_mtx(alignment_matrix, data = DMF5_antigen_table, # need to change data 
                           peptide_id_col = 'ContactPositions',  color_col = 'TCR',
                           alignment_threshold = 0,
                           vertex.label.degree = 0,  edge.weight = 0.2,
                           vertex.size= 5,  vertex.label.cex = 0.7)

## ---- fig.show='hold', fig.width=7, fig.height=6------------------------------
plot_network_distance_mtx(distance_mtx, data = DMF5_antigen_table,
                           peptide_id_col = 'ContactPositions',  color_col = 'TCR',
                           distance_threshold = 3,
                           vertex.label.degree = 0,  edge.weight = 0.2,
                           vertex.size= 5,  vertex.label.cex = 0.7, label_vertex=TRUE)

## ----citation-----------------------------------------------------------------
citation('diffSeqPatterns')

## ----sessionInfo, eval=TRUE---------------------------------------------------
sessionInfo()

