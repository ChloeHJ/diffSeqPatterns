#' This is selected peptide data from IEDB database
#'
#' @format A data frame with 1733 rows and 3 variables:
#' \describe{
#'   \item{Peptide}{Original 9aa peptide sequence}
#'   \item{Immunogenicity}{Positive or Negative denoting whether the peptide has been
#'   characterized to trigger immune response by T cell assays}
#'   \item{ContactPositions}{contact residues (positions 3-8) of the original peptide seqquence}
#' }
"iedbdata"

#' Positive peptide list
#'
#' This is list of peptide sequences annotated 'Positive' from T cell assays, analyzed for enrichment
#'
#' @format A list of 794 peptides
"pos_peptides"

#' Negative peptide list
#'
#' This is list of peptide sequences annotated 'Negative' from T cell assays, analyzed for depletion
#'
#' @format A list of 939 peptides
"neg_peptides"


