#' diffSeqPatterns: R package for analyzing differential peptide sequence patterns between two groups
#'
#' The diffSeqPatterns package contains functions to analyze and visualize
#' enrichment or depletion of sequence patterns between two groups,
#' such as amino acid usage, positional k-mer motifs, n-grams, and inter-sequence distances.
#'
#' It addresses the following needs:
#' \itemize{ \item Compute differential sequence patterns between two groups
#' \item Effectively visualize patterns;
#' \item Can be applied to any set of two peptide groups;
#' }
#'
#'
#' @section Section:
#' Common arguments to the functions.
#' @param peptides list of peptides
#' @param pos_peptides list of peptides analyzed for enrichment
#' @param neg_peptides list of peptides analyzed for depletion
#'
#' @import dplyr
#' @import caret
#' @import ggplot2
#' @import ngram
#' @import ggrepel
#' @import pheatmap
#' @import stringdist
#' @import ggseqlogo
#' @import PepTools
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices png pdf dev.off
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom scales muted
#' @importFrom data.table rbindlist
#' @importFrom stats predict
#' @importFrom stringr str_detect
#' @importFrom S4Vectors unique.Vector
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom igraph graph_from_adjacency_matrix E V plot.igraph
#' @importFrom reshape melt
#' @importFrom Biostrings pairwiseAlignment AAString AAStringSet
#'
#' @docType package
#' @name diffSeqPatterns
if(getRversion() >= "2.15.1")  utils::globalVariables(c( "pos_prop", "neg_prop", "pattern",
                                                         "count_neg",   "count_pos", "normCount_pos",
                                                         "normCount_neg", "ratio",
                                                         "ratioQual", 'counts', 'ngrams',
                                                         'label', 'position', 'aa',
                                                         'freq', 'prop', 'pos_freq',
                                                         'neg_freq', 'score', 'diff_pssm_mtx'
                                                         ))
NULL



#### To learn more about proBatch, start with the vignettes:
##### \code{browseVignettes(package = "diffSeqPatterns")}
