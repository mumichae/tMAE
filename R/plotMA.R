#' plotMA
#'
#' @description Creates an MA plot, ie, Fold Change (ALT/REF) vs Coverage colored 
#'    by significance (from p adjusted and allelic ratio)
#' @param data A data.frame containing the results table from \DESeq4MAE function
#' @param title The plot's title
#' @param padjCutoff The significance level
#' @param allelicRatioCutoff The minimum allelic ratio ALT/(ALT+REF) to be 
#'    considered signficant
#' @return A ggplot object containing the MA plot.
#' @export
#' @examples
#' file <- system.file("extdata", "demo_MAE_counts.tsv", package = "tMAE", mustWork = TRUE)
#' maeCounts <- read.table(file)
#' res <- run_deseq_all_mae(maeCounts)
#' plotMA(res)

plotMA <- function(data, title = NULL, padjCutoff = 0.05, 
                   allelicRatioCutoff = .8){
  stopifnot(c('altCount', 'refCount', 'altFreq', 'padj') %in% colnames(data))
  
  data <- as.data.table(data)
  data[, FC := (altCount + 1)/(refCount + 1)]
  data[, Significant := padj <= padjCutoff & 
         (altFreq > allelicRatioCutoff | altFreq < (1-allelicRatioCutoff))]
  ggplot(data, aes(totalCount, FC)) + geom_point(aes(col = Significant), size = .9) +
    theme_bw(base_size = 14) + scale_y_log10() + scale_x_log10() +
    scale_color_manual(values = c('gray61', 'firebrick')) +
    labs(x = 'RNA Coverage per heterozygous SNV', 
         y = 'Fold change of allelic counts\n(ALT+1)/(REF+1)', title = title)
}

# mt <- readRDS('/s/project/genetic_diagnosis/processed_results/mae/samples/65990-MUC1404_res.Rds')
# plotMA(rmae)

