#' Add allele frequencies from gnomAD
#'
#' @description Add allele frequency information from gnomAD.
#' @author Vicente Yepez
#' @param data A data.frame containing allelic counts.
#' @param gene_assembly either 'hg19' or 'hg38' indicating the genome assembly of the variants.
#' @param max_af_cutoff cutoff for a variant to be considered rare. Default is .001.
#' @param populations 
#' @return A data.table with the original contents plus columns containing allele frequencies from different gnomAD populations.
#' @export
#' @examples
#' file <- system.file("extdata", "demo_MAE_counts.tsv", package = "tMAE", mustWork = TRUE)
#' maeCounts <- read.table(file)
#' maeRes <- run_deseq_all_mae(maeCounts)
#' maeRes <- add_gnomAD_AF(maeCounts, gene_assembly = 'hg19')

add_gnomAD_AF <- function(data, gene_assembly = c('hg19', 'hg38'), max_af_cutoff = .001,
                      populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe')){
  mafdb <- switch('hg19',
                  hg19 = MafDb.gnomAD.r2.1.hs37d5,
                  hg38 = MafDb.gnomAD.r2.1.GRCh38
  )
  
  # Transform data into GRanges object
  pos <- data$position
  gr <- GRanges(seqnames = data$contig, ranges = IRanges(start = pos, end = pos), strand = '*')
  # Add score of all, African, American, East Asian and Non-Finnish European
  pt <- score(mafdb, gr, pop = populations)
  # Compute the MAX_AF
  pt$MAX_AF = apply(pt, 1, max, na.rm=TRUE)
  
  res <- cbind(data, pt) %>% as.data.table()
  
  # Replace Inf with NA
  res[is.infinite(MAX_AF), MAX_AF := NA]
  
  res[, rare := (MAX_AF <= max_af_cutoff | is.na(MAX_AF))]
  
  return(res)
}

# d2 <- add_gnomAD_AF(data)




