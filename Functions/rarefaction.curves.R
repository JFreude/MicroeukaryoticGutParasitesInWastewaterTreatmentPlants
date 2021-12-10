rarefaction.curves <- function(counts, step.size = 20){
  
  #' @title Rarefaction curves of microbial sequencing data
  #'
  #' @description This function forms a wrapper around the functions 'rarecurve' implemented in vegan, and returns their output in a format suitable for ggplot. 
  #' The rarefaction curves are calculated for each sample separately. Therefore, the total number of reads of one sample is split into increasing subsets. For each subset, the number of OTUs is counted. The size of the subsets is set by the argument \strong{\emph{'step.size'}}.
  #' 
  #' @section Required packages:
  #' vegan
  #'
  #' @param counts (\emph{data.frame, matrix}). A data.frame or matrix providing \strong{non-normalised} count data. Column names must be sample-IDs, row names taxa-IDs.
  #' @param step.size (\emph{integer}). One integer specifying the distance between the number of reads for which the number of OTUs is determined.
  #'
  #' @return data object
  #' @export
  #'
  #' @examples
  
  #############################################################################################
  #################################### Check prerequisites ####################################
  #----------------- Check if all needed functions and packages are available ----------------#
  if(!"vegan" %in% rownames(installed.packages())){
    stop("Please install the package 'vegan'.")
  }
  
  #--------------------------- Check if needed arguments are valid ---------------------------#
  # Check if 'counts' are provided
  if(missing(counts) | !(class(counts) %in% c("data.frame", "matrix"))){
    stop("Please provide a data.frame/matrix containing non-normalised count data with taxa as rows and samples as columns.")
  } else {
    # Transpose counts
    counts.transposed <- as.data.frame(t(counts))
  }

  # Check if the value of the argument 'step.size' is valid
  if(step.size %% 1 != 0 | class(step.size) != "numeric" | length(step.size) != 1){
    stop("Wrong value for the argument 'step.size'. Please provide one integer.")
  }
  
  #############################################################################################
  ################################ Rarefaction curves by sample ###############################

  # Calculate rarefaction curves per sample
  # Vegans rarecurves creates a plot, this one is save into an NULL pdf
  pdf(NULL)
  rarecurves <- vegan::rarecurve(counts.transposed, step = step.size)
  dev.off()
  
  # MAPPLY applies a function to all elements of a vector or a list, having the same position.
  # The function transforms data resulting from vegans function 'rarecurve' to a dataframe usable for ggplot.
  rarefaction.data <- mapply(FUN = function(rare.data, sample.ID){df <- as.data.frame(rare.data)
  colnames(df) <- "Number of OTUs"
  rownames(df) <- NULL
  `Number of reads` <- attr(rare.data, "Subsample")
  `Sample ID` <- sample.ID
  return(cbind(`Number of reads`,df,`Sample ID`))},
  rare.data = rarecurves, sample.ID = as.list(row.names(counts.transposed)), SIMPLIFY = FALSE)
  
  # Merge all datasets of the samples
  rarefaction.data <- do.call(rbind, rarefaction.data)
    
  return(rarefaction.data)
}
