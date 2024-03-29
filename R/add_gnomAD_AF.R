#' Add allele frequencies from gnomAD
#'
#' @description Add allele frequency information from gnomAD.
#' @author Vicente Yepez
#' @param data A data.frame containing allelic counts.
#' @param genome_assembly either 'hg19/hs37d5' or 'hg38/GRCh38' indicating the genome assembly of the variants.
#' @param max_af_cutoff cutoff for a variant to be considered rare. Default is .001.
#' @param populations The population to be annotated.
#' @param ... Used for backwards compatibility (gene_assembly -> genome_assembly)
#' @return A data.table with the original contents plus columns containing allele frequencies from different gnomAD populations.
#' @export
#' 
#' @examples
#' file <- system.file("extdata", "allelic_counts_HG00187.csv", package = "tMAE", mustWork = TRUE)
#' maeCounts <- fread(file)
#' maeRes <- DESeq4MAE(maeCounts)
#' maeRes <- add_gnomAD_AF(maeCounts, genome_assembly = 'hg19', pop="AF")
#' 
add_gnomAD_AF <- function(data, 
    genome_assembly = c('hg19', 'hs37d5', 'hg38', 'GRCh38'), 
    max_af_cutoff = .001,
    populations = c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_popmax'),
    ...){
  genome_assembly <- match.arg(genome_assembly)
  if("gene_assembly" %in% names(list(...))){
    warning("'gene_assembly' is deprecated. Please use 'genome_assembly' instead.")
    genome_assembly <- list(...)[['gene_assembly']]
  }
  if(genome_assembly %in% c('hg19', 'hs37d5')){
    if(!requireNamespace("MafDb.gnomAD.r2.1.hs37d5")){
      stop("Could not load gnomAD MafDb. Please install it.")
    }
    mafdb <- MafDb.gnomAD.r2.1.hs37d5::MafDb.gnomAD.r2.1.hs37d5
  } else if(genome_assembly %in% c('hg38', 'GRCh38')){
    if(!requireNamespace("MafDb.gnomAD.r2.1.GRCh38")){
      stop("Could not load gnomAD MafDb. Please install it.")
    }
    mafdb <- MafDb.gnomAD.r2.1.GRCh38::MafDb.gnomAD.r2.1.GRCh38
  } else {
    stop("Please provide a supported genome assembly version.")
  }
  
  if(!all(populations %in% populations(mafdb))){
    stop("Please provide only populations provided by gnomAD!")
  }
  
  # Transform data into GRanges object
  gr <- GRanges(seqnames = data$contig,
      ranges = IRanges(start=data$position, width=1), 
      strand = '*')
  
  # Add score of all, African, American, East Asian and Non-Finnish European
  pt <- score(mafdb, gr, pop = populations) %>% as.data.table()
  colnames(pt) <- populations
  res <- cbind(data, pt) %>% as.data.table()
  
  # Compute the MAX_AF (why do we change col names?)
  if(any(c("AF", "AF_popmax") %in% colnames(res))){
    res$MAX_AF <- apply(res[, ..populations], 1, max, na.rm = T)
    }
    
    # Replace Inf with NA
    res[is.infinite(MAX_AF), MAX_AF := NA]
    res[, rare := (MAX_AF <= max_af_cutoff | is.na(MAX_AF))]
  } else {
    warning("Could not find AF_popmax or AF column. Data is not filtered.")
  }
  
  return(res)
}

