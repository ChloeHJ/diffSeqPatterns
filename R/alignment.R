#' Compute pairwise alignment between peptides
#' @inheritParams diffSeqPatterns
#'
#' @param substitutionMatrix substitution matrix to compute pairwise alignment score,
#' examples include 'BLOSUM62', 'BLOSUM50' or can input matrix of interest e.g. \code{data('BLOSUM62')}
#' @param alignType type of alignment, options include "global", "local", "overlap",
#' "global-local", and "local-global"
#' @param gapOpening cost for opening a gap in the alignment.
#' @param gapExtension the incremental cost incurred along the length of the gap in the alignment.
#'
#' @return alignment matrix
#'
#' @examples
#' adjacency_matrix <- compute_pairwise_alignment(peptides = c(pos_peptides[1:50], neg_peptides[1:60]))
#' adjacency_matrix_positive <- compute_pairwise_alignment(peptides = pos_peptides[1:50])
#' adjacency_matrix_negative <- compute_pairwise_alignment(peptides = neg_peptides[1:60])
#'
#' @export
compute_pairwise_alignment <- function(peptides, substitutionMatrix = 'BLOSUM62', alignType = 'global',
                                       gapOpening=10, gapExtension=4){
  peptides <- unique(peptides)
  mtx <- matrix(NA, ncol =  length(peptides), nrow=length(peptides),
                dimnames=list(peptides,peptides))

  for(row in 1:nrow(mtx)){
    peptide_1 <- AAString(peptides[row])
    align <- lapply(peptides, function(elt_j) {
      pairwiseAlignment(peptide_1, AAString(elt_j), substitutionMatrix=substitutionMatrix,
                        type=alignType, gapOpening = gapOpening,
                        gapExtension = gapExtension, scoreOnly = TRUE)
    });
    mtx[row, ] <- unlist(align)
  }
  return(mtx)

}


#' Plot heatmap from pairwise alignment matrix
#' @param adjacency_matrix output from \code{compute_pairwise_alignment}
#' @param col_data dataframe that includes columns for color annotations and rownames of
#' peptides in adjacency matrix
#'
#' @return heatmap colored by distance between peptide sequences
#'
#' @examples
#' library(tidyverse)
#' adjacency_matrix <- compute_pairwise_alignment(peptides = c(pos_peptides[1:50], neg_peptides[1:60]))
#' adjacency_matrix_positive <- compute_pairwise_alignment(peptides = pos_peptides[1:50])
#'
#' col_data <- iedbdata %>% filter(ContactPositions %in% colnames(adjacency_matrix)) %>%
#'   select(ContactPositions, Immunogenicity) %>% distinct() %>%
#'   column_to_rownames(var =  'ContactPositions')
#' p1 <- plot_heatmap_alignment_mtx(adjacency_matrix, col_data)
#' p2 <- plot_heatmap_alignment_mtx(adjacency_matrix_positive, col_data)
#'
#' @export
plot_heatmap_alignment_mtx <- function(adjacency_matrix, col_data){
  make_symmetric_mtx <- function(mtx){
    b <- mtx
    b[upper.tri(b) ] <- t(b)[upper.tri(b )]

    return(b)
  }
  sym_mtx <- make_symmetric_mtx(adjacency_matrix)
  p <- pheatmap(sym_mtx, annotation_row = col_data, annotation_col = col_data)

  return(p)
}

#' Plot network graph from pairwise alignment matrix
#' @param adjacency_matrix output from \code{compute_pairwise_alignment}
#' @param data data to denote node colors - minimum columns are peptide and color_id
#' @param peptide_id_col column name for peptide_id e.g. Peptide
#' @param color_col column name for color annotation e.g. Immunogenicity
#' @param vertex.label.degree position of node label
#' @param edge.weight distance*edge.weight to represent thickness of edges in network plot
#' @param vertex.size node size
#' @param vertex.label.cex node font size
#' @param label_vertex TRUE/FALSE to print sequence for each node
#' @param alignment_threshold threshold to show e.g. alignment_threshold = 0 means
#' showing edges that have >= 0 alignment score
#'
#' @return Network graph showing distance between peptide sequences
#'
#' @examples
#' adjacency_matrix <- compute_pairwise_alignment(peptides = c(pos_peptides[1:50], neg_peptides[1:60]))
#' p1 <- plot_network_alignment_mtx(adjacency_matrix, data = NULL)
#' p2 <- plot_network_alignment_mtx(adjacency_matrix, data = iedbdata,
#'                            peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
#'                            vertex.label.degree = 0,  edge.weight = 0.2,
#'                            vertex.size= 5,  vertex.label.cex = 0.8, label_vertex=TRUE)
#' p3 <- plot_network_alignment_mtx(adjacency_matrix, data = iedbdata,
#'                            peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
#'                            edge.weight = 0.2,  vertex.size= 5,   label_vertex=FALSE)
#'
#' @export
plot_network_alignment_mtx <- function(adjacency_matrix, data = NULL,
                                       peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
                                       alignment_threshold = 0,
                                       vertex.label.degree = 0,  edge.weight = 0.2,
                                       vertex.size= 5,  vertex.label.cex = 0.8, label_vertex=FALSE){


  diag(adjacency_matrix) <- 0
  adjacency_matrix[adjacency_matrix<alignment_threshold] <- 0
  share.igraph <- graph_from_adjacency_matrix(adjacency_matrix, weighted=TRUE)

  # generate color annotation table
  if(!is.null(data)){
    data <- data[match(colnames(adjacency_matrix), data[[peptide_id_col]]),] %>%
      dplyr::select(!!c(peptide_id_col, color_col))

    getPalette = colorRampPalette(brewer.pal(6, "Dark2"))
    color_fac <- unique(data[[color_col]])
    color_dt <- data.frame(color = getPalette(length(color_fac)))
    rownames(color_dt) <- color_fac
    color_dt <- color_dt %>% tibble::rownames_to_column(var = color_col)
    color_table <- data %>% left_join(color_dt, by = c(color_col))

    igraph::V(share.igraph)$Node_color <- color_table$color
  }

  #igraph.layout <- layout_(share.igraph, with_dh(weight.edge.lengths = edge_density(share.igraph)/1000))
  igraph::E(share.igraph)$Edge_color <- c('grey')

  if(label_vertex == TRUE){
    p <- plot(share.igraph, vertex.label.color="black" , vertex.size=vertex.size,
              edge.width = (E(share.igraph)$weight)*edge.weight, edge.arrow.mode=0,
              vertex.color=V(share.igraph)$Node_color, edge.color = E(share.igraph)$Edge_color ,
              vertex.label.family= 'Helvetica', edge.curved	 = F, vertex.frame.color= 'White',
              vertex.label.cex = vertex.label.cex, #layout = igraph.layout,
              vertex.label.dist=1,  vertex.label.degree = vertex.label.degree)
  }


  if(label_vertex == FALSE){
    p <- plot(share.igraph, vertex.label.color="black" , vertex.size=vertex.size,
              edge.width = (E(share.igraph)$weight)*edge.weight, edge.arrow.mode=0,
              vertex.color=V(share.igraph)$Node_color, edge.color = E(share.igraph)$Edge_color ,
              vertex.label.family= 'Helvetica', edge.curved	 = F, vertex.frame.color= 'White',
              vertex.label=NA)
  }

  return(p)
}
