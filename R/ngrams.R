#' Compute enrichment or depletion of n-grams between two groups
#' @inheritParams diffSeqPatterns
#' @param ngram_lengths list of n to compute for n-grams (e.g. 2 for n-gram)
#'
#' @return data table summarizing number of n-grams appeared in positive and negative peptide lists
#' and ratio of counts normalized by the total number of positive and negative peptides respectively
#'
#' @examples  ngram_df <- compute_ngrams(pos_peptides, neg_peptides, ngram_lengths = c(2, 3, 4))
#'
#' @export
compute_ngrams <- function(pos_peptides, neg_peptides, ngram_lengths = c(2, 3, 4, 5)) {

  # split character string
  split_pos_peptides <- c()
  for(pep in pos_peptides){
    split_pos_peptides <-  c(split_pos_peptides, splitter(pep, split.char = TRUE ) )
  }
  split_neg_peptides <- c()
  for(pep in neg_peptides){
    split_neg_peptides <-  c(split_neg_peptides, splitter(pep, split.char = TRUE ) )
  }

  # compute n-grams and statistics
  ngram_table <- list()
  for(n in ngram_lengths){
    ng_pos <- get.phrasetable(ngram(split_pos_peptides, n=n)) %>% dplyr::rename(pos_freq = freq, pos_prop = prop)
    ng_neg <- get.phrasetable(ngram(split_neg_peptides, n=n)) %>% dplyr::rename(neg_freq = freq, neg_prop = prop)
    ngram_table[[n-1]] <- full_join(ng_pos,ng_neg, by = c('ngrams') )
  }; ngram_table <- rbindlist(ngram_table) %>% mutate(ratio = pos_prop/neg_prop)

  return(ngram_table)

}


#' Visualize enrichment or depletion score of n-grams in scatterplot
#' @param ngram_table output matrix from \code{compute_ngrams}
#' @param ratio_threshold ratio threshold to label n-grams
#' @param n_threshold for n-grams that are only present in either positive or negative list,
#' threshold to label n-grams, e.g. n_threshold = 4, only label motif when there are >= 4 peptides
#' containing the motif in either positive or negative peptide lists.
#' @return scatterplot
#'
#' @examples
#' ngram_df <- compute_ngrams(pos_peptides, neg_peptides, ngram_lengths = c(2, 3, 4))
#' p <- plot_point_ngrams(ngram_df)
#'
#' @export
plot_point_ngrams <- function(ngram_table,  ratio_threshold = 4,  n_threshold = 6){


  dat <- ngram_table %>% mutate(ngrams = gsub( " ", "", ngrams) ) %>%
    mutate(label = ifelse(is.infinite(ratio) & pos_freq >= n_threshold, ngrams,
                          ifelse(ratio == 0 & neg_freq >= n_threshold, ngrams,
                                 ifelse(!is.infinite(ratio) & ratio >= ratio_threshold, ngrams,
                                        ifelse(ratio != 0 & ratio <= 1/ratio_threshold, ngrams, ''))) ))

  p <- dat %>%
    ggplot(aes(x=pos_freq, y=neg_freq)) +
    geom_point() + theme_classic() +
    geom_text_repel(data = dat,
                    aes(label = label),
                    box.padding = 0.5, max.overlaps = Inf, #set max.overlap = Inf to show all labels
                    color = 'darkred',
                    show.legend = FALSE) +
    labs(x = 'Count in Positive', y = 'Count in Negative') +
    scale_x_continuous(trans='log10') +  scale_y_continuous(trans='log10')


  return(p)
}


#' Visualize top enriched and depleted n-grams in barplot
#' @param ngram_table output matrix from \code{compute_ngrams}
#' @param top_n number of top n-grams from positive and negative (respectively) to visualize
#' @param ratio_threshold ratio threshold to label n-grams
#' @param n_threshold for n-grams that are only present in either positive or negative list,
#' threshold to label n-grams, e.g. n_threshold = 4, only label motif when there are >= 4 peptides
#' containing the motif in either positive or negative peptide lists.
#' @return barplot
#'
#' @examples
#' ngram_df <- compute_ngrams(pos_peptides, neg_peptides, ngram_lengths = c(2, 3, 4))
#' p <- plot_bar_top_ngrams(ngram_df, top_n = 20)
#'
#' @export
plot_bar_top_ngrams <- function(ngram_table, top_n = 10,
                                ratio_threshold = 4,  n_threshold = 6){

  ngram_table <- ngram_table %>% mutate(ngrams = gsub( " ", "", ngrams) )
  dat <- ngram_table %>%
    mutate(ratioQual = ifelse(is.infinite(ratio) & pos_freq >= n_threshold, 'only_in_pos',
                              ifelse(ratio == 0 & neg_freq >= n_threshold, 'only_in_neg',
                                     ifelse(!is.infinite(ratio) & ratio >= ratio_threshold, 'high_in_pos',
                                            ifelse(ratio != 0 & ratio <= 1/ratio_threshold, 'high_in_neg', NA))) ))
  motifDF_select <- rbind(dat %>% filter(ratioQual == 'only_in_pos') %>% arrange(desc(pos_freq)) %>% head(n=top_n)  ,
                          dat %>% filter(ratioQual == 'high_in_pos') %>% arrange(desc(ratio)) %>% head(n=top_n)  ,
                          dat %>% filter(ratioQual == 'high_in_neg') %>% arrange(ratio) %>% head(n=top_n) ,
                          dat %>% filter(ratioQual == 'only_in_neg') %>% arrange(desc(neg_freq)) %>% head(n=top_n)) %>%
    dplyr::select(ngrams, pos_freq, neg_freq) %>% as.data.frame()

  order_motifs <- motifDF_select$ngrams
  data_plot <- motifDF_select %>% melt(id.vars = "ngrams"); colnames(data_plot) <- c('ngrams', 'class', 'counts')
  data_plot$ngrams <- factor(data_plot$ngrams, levels=order_motifs)
  data_plot <- data_plot %>%  mutate(class = ifelse(class == 'neg_freq', 'Negative', ifelse(class == 'pos_freq', 'Positive', NA))) %>%
    mutate(class = factor(class, levels = c('Positive', 'Negative')))

  g <- data_plot %>% ggplot(aes(x=ngrams, y=counts, fill=class)) +
    scale_fill_brewer(palette='Dark2') +
    geom_bar(stat="identity", width = 0.7)+theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5))

  return(g)

}

