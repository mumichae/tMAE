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
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicScores populations
#' @importFrom BiocGenerics score
#' @importFrom MafDb.gnomAD.r2.1.hs37d5 MafDb.gnomAD.r2.1.hs37d5
#' @importFrom MafDb.gnomAD.r2.1.GRCh38 MafDb.gnomAD.r2.1.GRCh38
#' @useDynLib tMAE
#' 
NULL

