beta.diversity<-function(counts, metadata, groups, subsample=TRUE, normalize=TRUE, 
                            wilcoxon.pvalues=TRUE, paired=FALSE, colored = TRUE, shaped=FALSE, return.indices=FALSE){
  
  #' @title Compare beta diversity across groups
  #'
  #' @description The function allows to compare the beta diversity across groups using the function vegdist implemented in vegan. 
  #' The beta diversities across groups is visualized in ggplot as box-plot. Optionally, a two-sided Wilcoxon test is conducted using the function stat_compare_means implemented in ggpubr, the significant results (p < 0.05) are shown in the plot.
  #'  
  #' @param counts (\emph{data.frame, matrix}). A data.frame or matrix providing count data. Column names must be sample-IDs, row names taxa-IDs.
  #' @param metadata (\emph{data.frame, matrix}). A data.frame or matrix providing metadata. Column names must be factors, row names must be sample-IDs.
  #' @param groups (\emph{character}). A character specifying one column name of the metadata, with either characters/factors or integers. The values of the selected columns are plotted on the x-axis of the resulting box plots.
  #' @param subsample (\emph{logical}). If \emph{TRUE}, groups are randomly subsampled to have the same sample number in each group, else data are not subsampled. By default, data are subsampled. 
  #' @param normalize (\emph{logical}). If \emph{TRUE}, count data are normalized row-wise (relative countss are calculated). By default, counts are normalized.
  #' @param property (\emph{"richness","evenness","alpha", "beta"}). Richness (Chao1), evenness (Sheldon), alpha (alpha-diversity with Shannon index) or beta (beta-diversity with Bray Curtis dissimilarities, only between-group pairs are used to calculate the beta diversity for all samples) diversities are calculated.
  #' @param wilcoxon.pvalues (\emph{logical}). If \emph{TRUE}, significant Wilcoxon p-values are displayed on the box plot using function stat_compare_means in R package ggpubr.
  #' @param paired (\emph{logical}). If \emph{TRUE} a paired, if \emph{FALSE} an unpaired wilcoxon test is conducted. By default an unpaired wilcoxon test is conducted.
  #' @param colored (\emph{logical/character}). If \emph{TRUE}, the samples are colored by the samples, if \emph{FALSE}, they are not colored. Alternatively, a \emph{character}, one column name of the \strong{metadata}, can be specified. By default the samples are colored by the samples.
  #' @param shaped (\emph{logical/character}). If \emph{TRUE}, the samples are shaped by the samples, if \emph{FALSE} they are not shaped. Alternatively, a \emph{character}, one column name of the \strong{metadata}, can be specified. By default the samples are not shaped.
  #' @param return.indices (\emph{logical}). If \emph{TRUE} a data frame with the diversity index is returned in addition to the plot (both stored in a list), if \emph{FALSE} only the plot is returned. By default return.indices = \emph{FALSE}.
  #'
  #' @section Packages the function bases on:
  #' 'ggplot2', 'ggpubr', 'reshape2', 'vegan'
  #'
  #' @section Required functions:
  #' 'calculate.beta.diversity', 'plot.beta.diversity'

  
  #############################################################################################
  #################################### Check prerequisites ####################################
  
  # Check if function 'calculate.beta.diversity' exists
  if(!file.exists(paste0(rprojroot::find_rstudio_root_file(), "/Functions/calculate.beta.diversity.R"))){
    stop("The function 'calculate.beta.diversity.R' could not be found in the directory 'Functions', please make sure that it is available.")
  } else {
    source(paste0(rprojroot::find_rstudio_root_file(), "/Functions/calculate.beta.diversity.R"))
  }
  
  # Check if function 'plot.beta.diversity' exists
  if(!file.exists(paste0(rprojroot::find_rstudio_root_file(), "/Functions/plot.beta.diversity.R"))){
    stop("The function 'plot.beta.diversity.R' could not be found in the directory 'Functions', please make sure that it is available.")
  } else {
    source(paste0(rprojroot::find_rstudio_root_file(), "/Functions/plot.beta.diversity.R"))
  }
  
  # Check if the value of the argument 'colored' is valid
  if(!colored %in% c(TRUE, FALSE, colnames(metadata))){
    stop("Wrong value for the argument 'colored' is provided. Please provide a logical (TRUE/FALSE) or a coloumn name of metadata.")
  }
  
  # Check if the value of the argument 'shaped' is valid
  if(!shaped %in% c(TRUE, FALSE, colnames(metadata))){
    stop("Wrong value for the argument 'shaped' is provided. Please provide a logical (TRUE/FALSE) or a coloumn name of metadata.")
  }
  
  # Check if the value of the argument 'return.indices' is valid
  if(class(return.indices) != "logical"){
    stop("Wrong value for the argument 'return.indices' is provided. Please provide a logical (TRUE/FALSE).")
  }
  
  #############################################################################################
  ################################# Calculate beta diversity ##################################

  # Set title and y axis label
  title = "Beta diversity"
  y.axis.label = "Bray-Curtis dissimilarity"
    
  # Calculate beta diversity
  diversity <- calculate.beta.diversity(counts=counts, metadata=metadata, groups=groups, subsample=subsample, normalize=normalize)
  
  # Prepare to return diversity table if desired
  if(return.indices){
    div <- cbind(Sample=rownames(diversity), diversity)
    rownames(div) <- NULL
    div <- cbind(Sample1 = gsub("\\..*$", "", div$Sample), Sample2 = gsub("^.*\\.", "", div$Sample),
                 div[,colnames(div) != "Sample"])
  }
  
  # Create column for colors
  if(!colored %in% c(T, F)){
    diversity <- setNames(cbind(diversity, metadata[match(gsub("\\..*$","",rownames(diversity)), rownames(metadata)), colored]),
                            c(colnames(diversity), colored))
  }
  
  # Create column for shapes
  if(!shaped %in% c(T, F)){
    # Check if shape and color have same names, if yes change it
    if(shaped == colored){
      colnames(metadata)[colnames(metadata) %in% shaped] <- paste(shaped, "1", sep = "_")
      shaped <- paste(shaped, "1", sep = "_")
    }
    
    diversity <- setNames(cbind(diversity, metadata[match(gsub("^.*\\.","",rownames(diversity)), rownames(metadata)), shaped]),
                          c(colnames(diversity), shaped))
  }
  
  # Plot diversity
  p <- plot.beta.diversity(plot.data=diversity, wilcoxon.pvalues=wilcoxon.pvalues, paired=paired, colored = colored, shaped = shaped) + ggplot2::theme_bw()
  
  
  # Add title and y axis label
  p <- p + ggplot2::labs(title=title, y=y.axis.label) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                   axis.title.x = ggplot2::element_text(size = 13),
                   axis.title.y = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_text(size = 13),
                   legend.text = ggplot2::element_text(size = 12))
  
  if(return.indices){
    results <- list(p, div)
    names(results) <- c("Plot", "Indices")
    return(results)
  } else {
    return(p)
  }
}
