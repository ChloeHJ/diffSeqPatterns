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


#' Positive peptide list from DMF5 selection
#'
#' This is list of peptide sequences at contact positions annotated 'Positive'
#' against DMF5 T cell screening, analyzed for enrichment.
#' Used for demonstrations in vignette and application note.
#'
#' @format A list of 55 peptides
"DMF5_pos_peptides"

#' Negative peptide list from DMF5 selection
#'
#' This is list of peptide sequences at contact positions annotated 'Negative'
#' from DMF5 T cell screening, analyzed for depletion.
#' Used for demonstrations in vignette and application note.
#'
#' @format A list of 200 peptides
"DMF5_neg_peptides"

#' Table of 55 DMF5 antigens and 200 background peptides annotated for their T cell specificity
#'
#' This is table used for color annotations. Used for demonstrations in vignette and application note.
#'
#' @format data table with peptides at their contact positions and their specificity against DMF5 T cells
"DMF5_antigen_table"
