# Converts a GRanges object to a Data Table
allelic_granges_to_dt <- function(data){
    # take only heterozygous mutation and with enough coverage
    data$GT <- as.character(data$GT)
    goodGT <- grepl('0[|/]1|1[|/]0', data$GT, perl=TRUE)
    # goodGT <- grepl('0[|/]1|1[|/]0', data$GT[,1], perl=TRUE)
    
    data$GQ <- as.integer(data$GQ)
    goodCov <- data$coverage > 0
    data <- data[goodGT & goodCov]
    
    data$nucl_piles = NULL
    data$qual_piles = NULL
    
    # get alt and ref counts
    data$chr <- seqnames(data)
    data$alt_cov <- as.integer(floor(data$coverage * data$alt_allele_freq))
    data$ref_cov <- as.integer(data$coverage - data$alt_cov)
    data$pos <- start(data)
    data$ALT <- as.character(unlist(data$ALT))
    dt <- as.data.table(data)
    dt[, c("start", "end", "width", "strand", "QUAL", "FILTER", "GQ", "A", "C", "G", "T", "seqnames") := NULL]
    return(dt)
}


# Performs the Negative Binomial Wald Test
deseq_for_allele_specific_expression <- function(data, minCoverage=10,
                                                     disp=0.05, independentFiltering=FALSE){
    # Obtain a data.table from imput
    if(any(class(data) == 'data.frame')){
        dt <- as.data.table(data)
        dt[, c('lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs') := NULL]
    } else if(class(data) == 'GRanges'){
        dt <- allelic_granges_to_dt(data)
    }
    
    dt <- dt[totalCount >= minCoverage]
    
    # create deseq object
    dds <- DESeqDataSetFromMatrix(
        as.matrix(dt[, .(altCount, refCount)]),
        DataFrame(condition=factor(c("altAllele", "refAllele"))),
        design = ~ condition
    )
    
    if(!is.null(data$hgncid))
        rownames(dds) <- data[, paste0(hgncid, "_", pos)]  # Doesn't always have the gene name
    
    mcols(rowRanges(dds)) <- dt
    
    # estimate the size factors and pseudo dispersion
    dds <- estimateSizeFactors(dds)
    # dds <- estimateDispersions(dds) # impossible to determine now with only 1 sample and 2 conditions
    
    # set dispersion by hand
    dispersions(dds) <- disp
    
    # run wald test
    dds <- nbinomWaldTest(dds)
    res <- results(dds, contrast = c("condition", "altAllele", "refAllele"), 
                   independentFiltering = independentFiltering)
    
    
    return(list(dt = dt, res = res))
}


# Get Allele Specific deseq results
get_allele_specific_deseq_results <- function(dds_res){
    
    # get needed info
    dt <- dds_res$dt
    
    # add pvalue and padj
    dt[, c("pvalue","padj", "log2FC") := list(dds_res$res$pvalue, dds_res$res$padj, dds_res$res$log2FoldChange)]
    
    # add frequency of the alternative allele
    dt[, altFreq := altCount / (altCount + refCount)]
    
    return(dt)
}

#' Run deseq test for MAE
#'
#' @description Uses a negative binomial test to determine if a variant is mono-allelically expressed.
#' @author Vicente Yepez, Christian Mertes
#' @param data A data.frame containing allelic counts.
#' @param minCoverage minimum total allelic count. Default is 10.
#' @param disp Gene dispersion for the NB test. Default is 0.05.
#' @param independentFiltering Parameter that affects the multiple testing. Default is FALSE.
#' @return Mono-allelic results table containing original counts plus p-value, p-adjusted and freqALT columns.
#' @export
#' @examples
#' file <- system.file("extdata", "demo_MAE_counts.tsv", package = "tMAE", mustWork = TRUE)
#' maeCounts <- read.table(file)
#' run_deseq_all_mae(maeCounts)

DESeq4MAE <- function(data, minCoverage = 10, disp = .05, independentFiltering = FALSE){
    
    pt <- deseq_for_allele_specific_expression(data, minCoverage=minCoverage, disp=disp, 
                                                   independentFiltering=independentFiltering)
    
    res <- get_allele_specific_deseq_results(pt)
    
    return(res)
}

# data <- fread('/data/ouga/home/ag_gagneur/yepez/workspace/RNAseq-ASHG19/Data/input_data/variants/chr21_allellic_counts.tsv')
# setnames(data, 'sample', 'MAE_ID')
# rmae <- DESeq4MAE(data)

