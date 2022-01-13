#' Compute differential position specific scoring matrix (PSSM) between two groups (applicable to peptides of same lengths)
#' @inheritParams diffSeqPatterns
#'
#' @return differnetial pssm matrix
#'
#' @examples diff_pssm_mtx <- compute_diff_pssm(pos_peptides, neg_peptides)
#'
#' @export
compute_diff_pssm <- function(pos_peptides, neg_peptides){

  all_peptides <- c(pos_peptides, neg_peptides)
  lengths <- nchar(all_peptides)
  if(length(unique(lengths)) != 1){stop('Not all sequences have equal lengths')}


  A_neg_PSSM <- Biostrings::consensusMatrix(AAStringSet(neg_peptides), as.prob = TRUE)
  A_pos_PSSM <- Biostrings::consensusMatrix(AAStringSet(pos_peptides), as.prob = TRUE)
  colnames(A_neg_PSSM) <- 1:ncol(A_neg_PSSM)
  colnames(A_pos_PSSM) <- 1:ncol(A_pos_PSSM)

  setdiff( rownames(A_pos_PSSM), rownames(A_neg_PSSM))
  A_pos_PSSM <- A_pos_PSSM[rownames(A_neg_PSSM), ]

  standardized_neg <- predict(preProcess(A_neg_PSSM, method=c("scale", 'center')), A_neg_PSSM)
  standardized_pos <- predict(preProcess(A_pos_PSSM, method=c("scale", 'center')), A_pos_PSSM)
  A_npos_pssm <- standardized_pos - standardized_neg

  colnames(A_npos_pssm) <- 1:ncol(A_npos_pssm)
  return(A_npos_pssm)
}


#' Plot differential PSSM into heatmap
#' @param diff_pssm_mtx output matrix from \code{compute_diff_pssm}
#'
#' @return PSSM plot
#'
#' @examples
#' diff_pssm_mtx <- compute_diff_pssm(pos_peptides, neg_peptides)
#' p <- plot_raster_diff_pssm(diff_pssm_mtx)
#'
#' @export
plot_raster_diff_pssm <- function(diff_pssm){
  longData<-melt(diff_pssm)
  colnames(longData) <- c('aa', 'position', 'score')
  pos <- 1:ncol(diff_pssm)

  p <- ggplot(longData, aes(x = position, y = aa)) +
    geom_raster(aes(fill=score)) +
    scale_fill_gradient2(low = muted('blue'), mid = 'white', high = muted('red')) +
    labs(x="Position", y="Amino acids") +
    scale_x_continuous(labels = as.character(pos), breaks = pos) +
    theme_minimal() + theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
                            axis.text.y=element_text(size=10),
                            plot.title=element_text(size=11, hjust = 0.5)) + theme_minimal()

  return(p)

}


#' Plot differential PSSM into sequence logos
#' @inheritParams diffSeqPatterns
#'
#' @return sequence logo plot
#'
#' @examples
#' diff_pssm_mtx <- compute_diff_pssm(pos_peptides, neg_peptides)
#' p <- plot_seqlogo_diff_pssm(pos_peptides, neg_peptides)
#'
#' @export
plot_seqlogo_diff_pssm <- function(pos_peptides, neg_peptides){
  pos_pssm <- pos_peptides %>%  pssm_freqs %>% pssm_bits %>% t
  neg_pssm <-  neg_peptides %>%  pssm_freqs %>% pssm_bits %>% t
  p <- (pos_pssm - neg_pssm) %>%  ggseqlogo(method="custom")
  return(p)
}









