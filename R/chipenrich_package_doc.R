#' chipenrich: Gene Set Enrichment For ChIP-seq Peak Data
#'
#' ChIP-Enrich performs gene set enrichment testing using peaks called from a
#' ChIP-seq experiment. The method empirically corrects for confounding factors
#' such as the length of genes, and the mappability of the sequence surrounding genes.
#'
#' @docType package
#' @name chipenrich
#'
#' @importFrom annotatr read_regions read_df
#' @import chipenrich.data
#' @import dplyr
#' @importClassesFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges findOverlaps
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @importFrom GenomeInfoDb seqinfo
#' @import grid
#' @importFrom IRanges IRanges
#' @import lattice
#' @import latticeExtra
#' @import methods
#' @import mgcv
#' @importFrom readr read_tsv
#' @import rms
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stats as.formula
#' @import stringr
#' @importFrom utils data
NULL
