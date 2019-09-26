#'
#' This is the import mapping for the tMAE package
#' 
#' @noRd
#' 
#' @name tMAE
#'
#' @import data.table
#' 
#' @importFrom BiocGenerics estimateSizeFactors plotDispEsts
#' 
#' @importFrom DESeq2 normalizationFactors normalizationFactors<-
#'          sizeFactors sizeFactors<- counts counts<-
#'          DESeqDataSetFromMatrix DESeqDataSet
#'          makeExampleDESeqDataSet show fpkm fpm
#'          estimateSizeFactorsForMatrix replaceOutliers
#'          dispersions dispersions<- nbinomWaldTest results
#' @useDynLib tMAE
#' 
NULL

