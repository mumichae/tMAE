test_that("End to end", {
    file <- system.file("extdata", "allelic_counts_HG00187.csv", 
                        package="tMAE", mustWork=TRUE)
    
    maeCounts <- fread(file)
    maeRes <- DESeq4MAE(maeCounts)
    
    expect_warning(add_gnomAD_AF(data=maeRes, gene_assembly='hg19', pop="AF"),
            "'gene_assembly' is depricated.")
})
