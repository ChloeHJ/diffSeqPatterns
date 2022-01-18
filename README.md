# diffSeqPatterns

diffSeqPatterns is an R package to comprehensively analyse differential sequence patterns between any two groups of peptides. The diffSeqPatterns provides functionality to analyse and visualize differential patterns in:

- Position-specific amino acid usage 
- N-grams: contiguous sequences of n-amino acids
- Position-specific k-mer motifs: contiguous or non-contiguous sequences of k amino acids
- Inter-sequence distance or alignment score: Pairwise distance (e.g. Hamming, levenshtein, Jaccard distance) or global/local alignment using substitution matrix of interest (e.g. BLOSUM62). 

Analysis and visualization tools in diffSeqPatterns can be applied to identify conserved patterns in any peptides, such as cancer neoepitopes, SARS-CoV-2 epitopes or antimicrobial peptides. 

## Installing

Install diffSeqPatterns from github:
```
library(devtools)
install_github("ChloeHJ/diffSeqPatterns", build_vignettes = TRUE)
```

Install diffSeqPatterns from CRAN:
```
library(diffSeqPatterns)
```


## Exploring the package

Browse the vignette:
```
browseVignettes('diffSeqPatterns')
```



