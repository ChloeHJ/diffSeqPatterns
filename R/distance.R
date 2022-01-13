#' Compute distance among peptides
#' @inheritParams diffSeqPatterns
#'
#' @param method Method for distance calculation, default is 'hamming'. Options include
#' "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", soundex".
#' See \code{'stringdist-metric'} for more details.
#'
#' @return distance matrix
#'
#' @examples
#' distance_matrix <- compute_pairwise_distance(peptides = c(pos_peptides, neg_peptides))
#'
#' @export
compute_pairwise_distance <- function(peptides, method = 'hamming'){

  peptides <- unique(peptides)
  mtx <- matrix(NA, ncol =  length(peptides), nrow=length(peptides),
                dimnames=list(peptides,peptides))

  for(row in 1:nrow(mtx)){
    peptide_1 <- peptides[row]
    align <- lapply(peptides, function(elt_j) {
      stringdist(peptide_1, elt_j, method = method)
    });
    mtx[row, ] <- unlist(align)
  }
  return(mtx)

}





#' Plot heatmap from distance adjacency matrix
#' @param distance_mtx output from \code{compute_pairwise_distance}
#' @param col_data dataframe that includes columns for color annotations and rownames of
#' peptides in adjacency matrix
#'
#' @return heatmap colored by distance between peptide sequences
#' @examples
#' distance_mtx <- compute_pairwise_distance(peptides = c(pos_peptides, neg_peptides))
#' col_data <- iedbdata %>% filter(ContactPositions %in% colnames(adjacency_matrix)) %>%
#'   select(ContactPositions, Immunogenicity) %>% distinct() %>% column_to_rownames(var =  'ContactPositions')
#' plot_heatmap_distance_mtx(distance_mtx, col_data)
#'
#' @export
plot_heatmap_distance_mtx <- function(distance_mtx, col_data){
  make_symmetric_mtx <- function(mtx){
    b <- mtx
    b[upper.tri(b) ] <- t(b)[upper.tri(b )]

    return(b)
  }
  sym_mtx <- make_symmetric_mtx(distance_mtx)
  sym_mtx <- sym_mtx*(-1)
  p <- pheatmap(sym_mtx, annotation_row = col_data, annotation_col = col_data)

  return(p)
}



#' Plot network graph from distance adjacency matrix
#' @param distance_mtx output from \code{compute_pairwise_distance}
#' @param data data to denote node colors - minimum columns are peptide and color_id
#' @param peptide_id_col column name for peptide_id e.g. Peptide
#' @param color_col column name for color annotation e.g. Immunogenicity
#' @param distance_threshold distance threshold to show e.g. distance_threshold = 3 means showing edges that have less or equal to 3 distance difference
#' @param vertex.label.degree position of node label
#' @param edge.weight distance*edge.weight to represent thickness of edges in network plot
#' @param vertex.size node size
#' @param vertex.label.cex node font size
#' @param label_vertex TRUE/FALSE to print sequence for each node
#'
#' @return Network graph showing distance between peptide sequences
#'
#' @examples
#' iedbdata <- read.csv(file = 'peptide_data.csv')
#' distance_mtx <- compute_pairwise_distance(peptides = c(pos_peptides[1:50], neg_peptides[1:60]))
#' p1 <- plot_network_distance_mtx(distance_mtx, data = iedbdata,
#'                           peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
#'                           distance_threshold = 3,
#'                           edge.weight = 0.3,  vertex.size= 5,   label_vertex=FALSE)
#' p2 <- plot_network_distance_mtx(distance_mtx, data = iedbdata,
#'                          peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
#'                           distance_threshold = 3,
#'                           edge.weight = 0.3,  vertex.size= 5,   label_vertex=TRUE)
#' p3 <- plot_network_distance_mtx(distance_mtx, data = iedbdata,
#'                           peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
#'                           distance_threshold = 4,
#'                           edge.weight = 0.3,  vertex.size= 5,   label_vertex=FALSE)
#'
#' @export
plot_network_distance_mtx <- function(distance_mtx, data = NULL,
                                      peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
                                      distance_threshold = 3,
                                      vertex.label.degree = 0,  edge.weight = 0.2,
                                      vertex.size= 5,  vertex.label.cex = 0.8, label_vertex=FALSE){


  distance_mtx <- 1/distance_mtx
  distance_mtx[distance_mtx< (1/distance_threshold)] <- 0
  distance_mtx[is.infinite(distance_mtx)] <- 0

  share.igraph <- graph_from_adjacency_matrix(distance_mtx, weighted=TRUE)

  # generate color annotation table
  if(!is.null(data)){
    data <- data[match(colnames(distance_mtx), data[[peptide_id_col]]),] %>%
      dplyr::select(!!c(peptide_id_col, color_col))

    getPalette = colorRampPalette(brewer.pal(6, "Dark2"))
    color_fac <- unique(data[[color_col]])
    color_dt <- data.frame(color = getPalette(length(color_fac)))
    rownames(color_dt) <- color_fac
    color_dt <- color_dt %>% rownames_to_column(var = color_col)
    color_table <- data %>% left_join(color_dt, by = c(color_col))

    V(share.igraph)$Node_color <- color_table$color
  }

  #igraph.layout <- layout_(share.igraph, with_dh(weight.edge.lengths = edge_density(share.igraph)/1000))
  E(share.igraph)$Edge_color <- c('grey')

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

