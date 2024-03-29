% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance.R
\name{plot_network_distance_mtx}
\alias{plot_network_distance_mtx}
\title{Plot network graph from distance adjacency matrix}
\usage{
plot_network_distance_mtx(
  distance_mtx,
  data = NULL,
  peptide_id_col = "ContactPositions",
  color_col = "Immunogenicity",
  distance_threshold = 3,
  vertex.label.degree = 0,
  edge.weight = 0.2,
  vertex.size = 5,
  vertex.label.cex = 0.8,
  label_vertex = FALSE
)
}
\arguments{
\item{distance_mtx}{output from \code{compute_pairwise_distance}}

\item{data}{data to denote node colors - minimum columns are peptide and color_id}

\item{peptide_id_col}{column name for peptide_id e.g. Peptide}

\item{color_col}{column name for color annotation e.g. Immunogenicity}

\item{distance_threshold}{distance threshold to show e.g. distance_threshold = 3 means
showing edges that have less or equal to 3 distance difference}

\item{vertex.label.degree}{position of node label}

\item{edge.weight}{distance*edge.weight to represent thickness of edges in network plot}

\item{vertex.size}{node size}

\item{vertex.label.cex}{node font size}

\item{label_vertex}{TRUE/FALSE to print sequence for each node}
}
\value{
Network graph showing distance between peptide sequences
}
\description{
Plot network graph from distance adjacency matrix
}
\examples{
distance_mtx <- compute_pairwise_distance(peptides = c(pos_peptides[1:50], neg_peptides[1:60]))
p1 <- plot_network_distance_mtx(distance_mtx, data = iedbdata,
                          peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
                          distance_threshold = 3,
                          edge.weight = 0.3,  vertex.size= 5,   label_vertex=FALSE)
p2 <- plot_network_distance_mtx(distance_mtx, data = iedbdata,
                         peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
                          distance_threshold = 3,
                          edge.weight = 0.3,  vertex.size= 5,   label_vertex=TRUE)
p3 <- plot_network_distance_mtx(distance_mtx, data = iedbdata,
                          peptide_id_col = 'ContactPositions',  color_col = 'Immunogenicity',
                          distance_threshold = 4,
                          edge.weight = 0.3,  vertex.size= 5,   label_vertex=FALSE)

}
