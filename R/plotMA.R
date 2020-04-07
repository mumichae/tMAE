#' plotMA4MAE
#'
#' @description Creates an MA plot, ie, Fold Change (ALT/REF) vs Coverage colored 
#'    by significance (from p adjusted and allelic ratio) or significance and 
#'    minor allele frequency (if rare_column is provided).
#' @author Vicente Yepez
#' @param data A data.frame containing the results table from \code{DESeq4MAE} function
#' @param title The plot's title
#' @param padjCutoff The significance level
#' @param allelicRatioCutoff The minimum allelic ratio ALT/(ALT+REF) to be 
#'    considered signficant
#' @param rare_column The name of the column that indicates if a variant is rare or not.
#'    Default is \code{null} which means it won't be plotted.
#' @return A ggplot object containing the MA plot.
#' @export
#' @examples
#' # file <- system.file("extdata", "demo_MAE_counts.tsv",
#'  # package = "tMAE", mustWork = TRUE)
#' # maeCounts <- read.table(file)
#' # res <- run_deseq_all_mae(maeCounts)
#' # plotMA4MAE(res)

plotMA4MAE <- function(data, title = NULL, padjCutoff = 0.05, 
                   allelicRatioCutoff = .8, rare_column = NULL){
  stopifnot(c('altCount', 'refCount', 'altRatio', 'padj') %in% colnames(data))
  
  data <- as.data.table(data)
  data[, FC := (altCount + 1)/(refCount + 1)]
  data[, Significant := padj <= padjCutoff & 
         (altRatio > allelicRatioCutoff | altRatio < (1-allelicRatioCutoff))]
  
  # Make the sketch of the plot
  g <- ggplot(data, aes(totalCount, FC)) + 
    theme_bw(base_size = 14) + scale_y_log10() + scale_x_log10() +
    labs(x = 'RNA Coverage per heterozygous SNV', 
         y = 'Fold change of allelic counts\n(ALT+1)/(REF+1)', title = title)
  
  # If rare_column provided, combine it with Significant column and plot
  if(!is.null(rare_column)){
    stopifnot(rare_column %in% colnames(data))
    data[, class := 'NS']
    data[Significant == T & get(rare_column) == T, class := 'Significant\n& Rare']
    data[Significant == T & get(rare_column) == F, class := 'Significant']
    g <- g + geom_point(aes(col = class), size = .9) +
      scale_color_manual(values = c('gray61', 'chocolate1', 'firebrick')) +
      theme(legend.title = element_blank())
  } else{
    g <- g + geom_point(aes(col = Significant), size = .9) +
      scale_color_manual(values = c('gray61', 'firebrick'))
  }
  
  return(g)
}

# mt <- readRDS('/s/project/genetic_diagnosis/processed_results/mae/samples/65990-MUC1404_res.Rds')
# plotMA(rmae, rare_column = 'rare')

#' plotAllelicCounts
#'
#' @description Creates an allelic counts plot, ie, counts of the alternative
#'    vs counts of the reference on the log10 scale. If significant and minor
#'    allele frequency columns present, dots will be highlighted accordingly.
#' @author Vicente Yepez
#' @param data A data.frame containing the counts or results table from \code{DESeq4MAE} function
#' @param title The plot's title
#' @param padjCutoff The significance level
#' @param allelicRatioCutoff The minimum allelic ratio ALT/(ALT+REF) to be 
#'    considered signficant
#' @param rare_column The name of the column that indicates if a variant is rare or not.
#'    Default is \code{null} which means it won't be plotted.
#' @return A ggplot object containing the MA plot.
#' @export
#' @examples
#' # file <- system.file("extdata", "demo_MAE_counts.tsv",
#'  # package = "tMAE", mustWork = TRUE)
#' # maeCounts <- read.table(file)
#' # res <- run_deseq_all_mae(maeCounts)
#' # plotAllelicCounts(res)

plotAllelicCounts <- function(data, title = NULL, padjCutoff = 0.05, 
                       allelicRatioCutoff = .8, rare_column = NULL){
  stopifnot(c('altCount', 'refCount') %in% colnames(data))
  
  data <- as.data.table(data)
  
  # Make the sketch of the plot
  g <- ggplot(data, aes(refCount, altCount)) + 
    theme_bw(base_size = 14) + scale_y_log10() + scale_x_log10() +
    labs(x = 'REF count + 1', y = 'ALT count + 1', title = title) + 
    geom_abline(slope = 1, intercept = 0)
  
  # If significance column is provided, add it to the plot
  if(!is.null(data$padj)){
    data[, Significant := padj <= padjCutoff & 
           (altRatio > allelicRatioCutoff | altRatio < (1-allelicRatioCutoff))]
  }
  
  # If rare_column provided, combine it with Significant column and plot
  if(!is.null(rare_column)){
    stopifnot(rare_column %in% colnames(data))
    data[, class := 'NS']
    data[Significant == T & get(rare_column) == T, class := 'Significant\n& Rare']
    data[Significant == T & get(rare_column) == F, class := 'Significant']
    g <- g + geom_point(aes(col = class), size = .9) +
      scale_color_manual(values = c('gray61', 'chocolate1', 'firebrick')) +
      theme(legend.title = element_blank())
  } else{
    g <- g + geom_point(aes(col = Significant), size = .9) +
      scale_color_manual(values = c('gray61', 'firebrick')) 
  }
  return(g)
}

# plotAllelicCounts(resMAE)
# plotAllelicCounts(resMAE, rare_column = 'rare')
