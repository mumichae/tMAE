test_that("End to end", {
    file <- system.file("extdata", "allelic_counts_HG00187.csv", 
            package="tMAE", mustWork=TRUE)
    
    maeCounts <- fread(file)
    maeRes <- DESeq4MAE(maeCounts)
    
    if(!requireNamespace("MafDb.gnomAD.r2.1.GRCh38", quietly=TRUE)){
        expect_error(add_gnomAD_AF(data=maeCounts, genome_assembly='hg38'), 
                "Could not load gnomAD MafDb")
        expect_error(add_gnomAD_AF(data=maeCounts, genome_assembly='GRCh38'),
                "Could not load gnomAD MafDb")
    } else {
        res <- add_gnomAD_AF(data=maeCounts, genome_assembly='hg38')
        res <- add_gnomAD_AF(data=maeCounts, genome_assembly='GRCh38')
    }
    
    if(!requireNamespace("MafDb.gnomAD.r2.1.hs37d5", quietly=TRUE)){
        if (!requireNamespace("BiocManager", quietly=TRUE))
            install.packages("BiocManager")
        BiocManager::install("MafDb.gnomAD.r2.1.hs37d5")
    }
    
    res <- add_gnomAD_AF(data=maeCounts, genome_assembly='hg19', pop=c("AF"))
        

    expect_is(plotMA4MAE(res), "ggplot")
    expect_is(plotAllelicCounts(res), "ggplot")
})
