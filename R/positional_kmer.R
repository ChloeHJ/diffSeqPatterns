#' Compute enrichment or depletion of positional k-mers between two groups (applicable to peptides of same lengths)
#' @inheritParams diffSeqPatterns
#'
#' @return data table summarizing number of positional k-mer motifs appeared in postiive and negative peptide lists
#' and ratio of counts normalized by the total number of positive and negative peptides respectively
#'
#' @examples kmer_df <- compute_positional_kmers(pos_peptides, neg_peptides)
#'
#' @export
compute_positional_kmers <- function(pos_peptides, neg_peptides){

  generate_motifs <- function(kmer_seqVector, k){
    require(S4Vectors)

    #1-mer
    # aa <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
    # k1mer <- NULL
    # for(pos in 1:k){
    #   for(i in 1:length(aa)){
    #     plain <- paste(c(rep('.', k) ), sep="",collapse="")
    #     substring(plain, pos) <- aa[i]
    #     k1mer <- c(k1mer, plain)
    #   }
    # }

    kmer <- list()
    for(kmer_gap in 1:(k-1)){
      kmer[[kmer_gap]] <- list()

      for(gap in kmer_gap:(k-1)){  #designate gap depending on kmer
        kmer[[kmer_gap]][[gap]] <- list()

        for(pos in 1:(k-gap)){ #scanning through position
          kmer[[kmer_gap]][[gap]][[pos]] <- unique(substr(kmer_seqVector, start=pos, stop=pos + gap));

          if(gap > kmer_gap){ # within position, where gap needs to be introduced within motif
            all_gapped_strings <- c()
            for(pos_to_gap in 2:gap){ # gaps in all possible conditions
              strings <-  kmer[[kmer_gap]][[gap]][[pos]]
              substring(strings, pos_to_gap, pos_to_gap-2+gap-(kmer_gap -1)) <- paste(c(rep('.', gap-(kmer_gap-1)-1)), sep="", collapse="")
              all_gapped_strings <- c(all_gapped_strings, strings)
            };kmer[[kmer_gap]][[gap]][[pos]] <- all_gapped_strings
          }

          kmer[[kmer_gap]][[gap]][[pos]]  <- paste0(paste(c(rep('.', pos-1 )), sep="",collapse=""), kmer[[kmer_gap]][[gap]][[pos]], paste(c(rep('.', k- (pos+gap) )), sep="",collapse=""))
        }
      }
    }
    kmer <- unlist(lapply(kmer, unique.Vector))

    # 5-mer
    nmer <- unique(kmer_seqVector)

    kmer_pattern <- as.character(c(kmer, nmer))
    return(kmer_pattern)

  }
  count_motifs <- function(motifs, pos_peptides, neg_peptides, enrichment_threshold, sum_threshold){

    motifDF <- data.frame(pattern = motifs, count_neg = NA, count_pos = NA) %>% mutate(pattern = as.character(pattern))
    for(row in 1:nrow(motifDF)){
      motifDF[row, "count_neg"] <- sum(str_detect(neg_peptides, motifDF[row,"pattern"]))
      motifDF[row, "count_pos"] <- sum(str_detect(pos_peptides, motifDF[row,"pattern"]))
    }

    motifDF <- motifDF %>%
      mutate(normCount_neg = count_neg/length(neg_peptides)) %>%
      mutate(normCount_pos = count_pos/length(pos_peptides)) %>%
      mutate(ratio = normCount_pos/normCount_neg) %>% arrange(-normCount_pos) %>% arrange(-ratio)


    return(motifDF)

  }

  all_peptides <- c(pos_peptides, neg_peptides)
  lengths <- nchar(all_peptides)
  if(length(unique(lengths)) != 1){stop('Not all sequences have equal lengths')}
  pep_length <- unique(lengths)

  motifs <- generate_motifs(all_peptides %>% unique,  k = pep_length) %>% unique
  motif_df <- count_motifs(motifs,  pos_peptides,  neg_peptides)
  return(motif_df)

  }


#' Visualize enrichment or depletion score of k-mer motifs in scatterplot
#' @param positional_kmer_data output matrix from \code{compute_positional_kmers}
#' @param ratio_threshold ratio threshold to label k-mer motifs
#' @param n_threshold for k-mer motifs that are only present in either positive or negative list,
#' threshold to label k-mer motifs, e.g. n_threshold = 4, only label motif when there are >= 4 peptides
#' containing the motif in either positive or negative peptide lists.
#'
#' @return scatterplot
#'
#' @examples
#' positional_kmer_data <- compute_positional_kmers(pos_peptides, neg_peptides)
#' p <- plot_point_positional_kmers(positional_kmer_data)
#'
#'
#' @export
plot_point_positional_kmers <- function(positional_kmer_data,
                                        ratio_threshold = 4,  n_threshold = 6){

   dat <- positional_kmer_data %>%
    mutate(label = ifelse(is.infinite(ratio) & count_pos >= n_threshold, pattern,
                          ifelse(ratio == 0 & count_neg >= n_threshold, pattern,
                                 ifelse(!is.infinite(ratio) & ratio >= ratio_threshold, pattern,
                                        ifelse(ratio != 0 & ratio <= 1/ratio_threshold, pattern, ''))) ))
   p <- dat %>%
    ggplot(aes(x=count_pos, y=count_neg)) +
    geom_point() + theme_classic() +
    geom_text_repel(data = dat,
                    aes(label = label),
                    box.padding = 0.5, max.overlaps = Inf, #set max.overlap = Inf to show all labels
                    color = 'darkred',
                    show.legend = FALSE) +
     labs(x = 'Count in Positive', y = 'Count in Negative')


  return(p)
}


#' Visualize top enriched and depleted k-mer motifs in barplot
#' @param positional_kmer_data output matrix from \code{compute_positional_kmers}
#' @param top_n number of top k-mer motifs from positive and negative (respectively) to visualize
#' @param ratio_threshold ratio threshold to label k-mer motifs
#' @param n_threshold for k-mer motifs that are only present in either positive or negative list,
#' threshold to label k-mer motifs, e.g. n_threshold = 4, only label motif when there are >= 4 peptides
#' containing the motif in either positive or negative peptide lists.
#'
#' @return barplot
#'
#' @examples
#' positional_kmer_data <- compute_positional_kmers(pos_peptides, neg_peptides)
#' plot_bar_posiitonal_kmers(positional_kmer_data, top_n = 20)
#'
#' @export
plot_bar_posiitonal_kmers <- function(positional_kmer_data, top_n = 10,
                                      ratio_threshold = 4,  n_threshold = 6){

  dat <- positional_kmer_data %>%
    mutate(ratioQual = ifelse(is.infinite(ratio) & count_pos >= n_threshold, 'only_in_pos',
                          ifelse(ratio == 0 & count_neg >= n_threshold, 'only_in_neg',
                                 ifelse(!is.infinite(ratio) & ratio >= ratio_threshold, 'high_in_pos',
                                        ifelse(ratio != 0 & ratio <= 1/ratio_threshold, 'high_in_neg', NA))) ))
  motifDF_select <- rbind(dat %>% filter(ratioQual == 'only_in_pos') %>% arrange(desc(count_pos)) %>% head(n=top_n)  ,
                          dat %>% filter(ratioQual == 'high_in_pos') %>% arrange(desc(ratio)) %>% head(n=top_n)  ,
                          dat %>% filter(ratioQual == 'high_in_neg') %>% arrange(ratio) %>% head(n=top_n) ,
                          dat %>% filter(ratioQual == 'only_in_neg') %>% arrange(desc(count_neg)) %>% head(n=top_n)) %>%
    dplyr::select(pattern, count_neg, count_pos) %>% as.data.frame()

  order_motifs <- motifDF_select$pattern
  data_plot <- motifDF_select %>% melt(id.vars = "pattern"); colnames(data_plot) <- c('pattern', 'class', 'counts')
  data_plot$pattern <- factor(data_plot$pattern, levels=order_motifs)
  data_plot <- data_plot %>%  mutate(class = ifelse(class == 'count_neg', 'Negative', ifelse(class == 'count_pos', 'Positive', NA))) %>%
    mutate(class = factor(class, levels = c('Positive', 'Negative')))

  g <- data_plot %>% ggplot(aes(x=pattern, y=counts, fill=class)) +
    scale_fill_brewer(palette='Dark2') +
    geom_bar(stat="identity", width = 0.7)+theme_bw() +
    theme(axis.text.x = element_text(angle = 90,  hjust=1),
          plot.title = element_text(hjust = 0.5))

  print(g)

}
