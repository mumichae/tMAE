#' Convert GRanges object to a Data Table
#'
#' This function allows you to convert a GRanges object to a data table.
#' @param data GRanges object
#' @return data.table object
#' @examples
#' allelic_granges_to_dt()

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


#' Run Deseq for allele Specific Expression
#'
#' This function allows you to ...
#' @param data data.table object.
#' @param min_cov Default is 10.
#' @param disp Default is 0.05
#' @param indepFilter Default is False
#' @return A list with the results.
#' @examples
#' run_deseq_for_allele_specific_expression()

run_deseq_for_allele_specific_expression <- function(data, min_cov=10,
                    disp=0.05, indepFilter=FALSE){

    require(DESeq2)

    if(any(class(data) == 'data.table')){
        dt <- copy(data)
        dt[, c('lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs') := NULL]
        setnames(dt, old = c("contig", "refAllele", "altAllele", "refCount", "altCount", "totalCount", "position"),
                     new = c("chr", "REF", "ALT", "ref_cov", "alt_cov", "coverage", "pos"))
    } else if(class(data) == 'GRanges'){
        dt <- allelic_granges_to_dt(data)
    }

    dt <- dt[coverage >= min_cov]

    # create deseq object
    dds <- DESeqDataSetFromMatrix(
        as.matrix(dt[, .(alt_cov, ref_cov)]),
        DataFrame(condition=factor(c("alt", "ref"))),
        design = ~ condition
    )

    if(!is.null(data$hgncid))
        rownames(dds) <- data[, paste0(hgncid, "_", pos)]  # Doesn't always have the gene name

    mcols(rowRanges(dds)) <- dt

    # estimate the size factor and pseudo dispersion
    dds <- estimateSizeFactors(dds)
    # dds <- estimateDispersions(dds) # impossible to determine now with only 1 sample and 2 conditions

    # set dispersion by hand
    dispersions(dds) <- disp

    # run wald test
    dds <- nbinomWaldTest(dds)
    res <- results(dds, contrast = c("condition", "alt", "ref"), independentFiltering = indepFilter)


    return(list(dds = dds, res = res))
}


#' Get Allele Specific deseq results
#'
#' This function allows you to ...
#' @param dds_res object
#' @param min_fc Default is log2(3/1)
#' @param padj_threshold Default is 0.1
#' @return data.table object
#' @examples
#' get_allele_specific_deseq_results()

get_allele_specific_deseq_results <- function(dds_res, min_fc=log2(3/1), padj_threshold = .1){

    # get needed info
    dds <- dds_res$dds
    mc <- mcols(dds)
    mc$REF <- as.character(mc$REF)
    mc$ALT <- as.character(unlist(mc$ALT))
    mc$GT = '0/1'
    dt <- as.data.table(mc[, c("sample", "chr", "pos", "REF", "ALT", "GT", "alt_cov", "ref_cov")])

    # add pvalue and padj
    dt[, c("pvalue","padj", "log2FC") := list(dds_res$res$pvalue, dds_res$res$padj, dds_res$res$log2FoldChange)]

    # prediction
    dt[, as_gt := "0/1"]
    dt[padj <= padj_threshold & abs(log2FC) >= min_fc, as_gt := ifelse(log2FC < 0, "0/0", "1/1")]

    # add additional features
    dt[, alt_freq := alt_cov / (alt_cov + ref_cov)]

    return(dt)
}

#' Run deseq test for MAE
#'
#' This function allows you to ...
#' @param data GRanges object.
#' @param min_cov Default is 10
#' @param disp Default is 0.05
#' @param indepFilter Default is False
#' @param min_fc Default is log2(3/1)
#' @param padj_threshold Default is 0.01
#' @return Results.
#' @export
#' @examples
#' run_deseq_all_mae()

run_deseq_all_mae <- function(data, min_cov = 10, disp = .05, indepFilter = FALSE, min_fc=log2(3/1), padj_threshold = .1){

    pt <- run_deseq_for_allele_specific_expression(data, min_cov=min_cov, disp=disp, indepFilter=indepFilter)

    res <- get_allele_specific_deseq_results(pt, min_fc = min_fc, padj_threshold = padj_threshold)

    return(res)
}

