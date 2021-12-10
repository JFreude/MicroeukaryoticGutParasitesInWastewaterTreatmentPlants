calculate.NMDS <- function(counts, dis = "bray", dimensions = 2, maxit=1000, trymax=100, ...){
  
  #' @title Non-metric Multidimensional Scaling for microbial sequencing data
  #'
  #' @description This function forms a wrapper around the function 'metaMDS' implemented in vegan, returning results usable with ggplot. 
  #' The function returns a list containing 'Eigenvectors', 'EigenvectorsSpecies', 'StressValue' and the results of the NMDS 'Ordination'. 
  #' 'Eigenvectors', 'EigenvectorsSpecies' and 'Stress value' can be used to visualize data with ggplot. 'Ordination' can be used to fit environmental vectors or factors onto the ordination using vegans 'envfit' function.
  #'
  #' @section Required packages:
  #' 'vegan'
  #'
  #' @param counts (\emph{data.frame, matrix}). A data.frame or matrix providing count data. Column names must be sample-IDs, row names taxa-IDs. Data must not contain any NAs. 
  #' @param dis (\emph{"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis", "chisq", "chord"}). Dissimilarity or distance indices supported by vegan, that can be used to calculate NMDS.
  #' @param dimensions (\emph{integer}). The number of dimensions used for the NMDS analysis. The argument must be given if the argument \strong{\emph{'ordination.methods'}} is set to \strong{"NMDS"}.
  #' @param maxit (\emph{integer}). Maximum number of iterations in the single NMDS run (see \code{\link[vegan]{metaMDS}}).
  #' @param trymax (\emph{integer}). Maximum numbers of random starts in search of stable solution (see \code{\link[vegan]{metaMDS}}).
  #' @param ... arguments passed to \code{\link[vegan]{metaMDS}}.
  #'
  #' @return
  #' @export
  #'
  #' @examples

  #############################################################################################
  #################################### Check prerequisites ####################################
  #--------------------------- Check if needed arguments are valid ---------------------------#
  
  # Check if 'counts' are provided
  if(missing(counts) | !(class(counts) %in% c("data.frame", "matrix"))){
    stop("Please provide count data formatted as data.frame or matrix")
  } else if(!all(unlist(lapply(counts, is.numeric)))){
    stop("Please check your count data, some of the entries are not numeric")
  } else if(class(counts) == "matrix"){
    counts <- as.data.frame(counts)
  }

  # Check if the value of the argument 'dis' is valid
  if(!dis %in% c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis", "chisq", "chord")){
    stop("Wrong value for the argument 'dis' is provided. Please choose between 'manhattan', 'euclidean', 'canberra', 'bray', 'kulczynski', 'jaccard', 'gower', 'altGower', 'morisita', 'horn', 'mountford', 'raup' , 'binomial', 'chao', 'cao', 'mahalanobis', 'chisq' and 'chord'.")
  }

  # Check if the value of the argument 'dimensions' is valid
  if(!is.numeric(dimensions)){# check if numeric
    stop("Please provide an integer for the argument 'dimensions'.")
  } else if(dimensions%%1 != 0){# check if integer
    stop("Please provide an integer for the argument 'dimensions'.")
  }
  
  # Check if the value of the argument 'maxit' is valid
  if(!is.numeric(maxit)){# check if numeric
    stop("Please provide an integer for the argument 'maxit'.")
  } else if(maxit%%1 != 0){# check if integer
    stop("Please provide an integer for the argument 'maxit'.")
  }
  
  # Check if the value of the argument 'trymax' is valid
  if(!is.numeric(trymax)){# check if numeric
    stop("Please provide an integer for the argument 'trymax'.")
  } else if(trymax%%1 != 0){# check if integer
    stop("Please provide an integer for the argument 'trymax'.")
  }
  
  #############################################################################################
  ####################################### Carry out NMDS ###################################### 
    
  # Delete samples with no entries
  if(any(colSums(counts)==0)){
    warning(paste0("No taxa are present in the sample(s) ", paste(names(which(colSums(counts)==0)), collapse = " & "),
                   ", thus they were excluded from the analysis"))
    counts <- counts[,colnames(counts)[colSums(counts)!=0]]
  }
  
  # Produces replicable results
  set.seed(12345)
  
  # Carry out NMDS
  NMDS <- vegan::metaMDS(data.frame(t(counts)), dist=dis, k=dimensions, maxit=maxit, trymax=trymax, ...)

  # Extract eigenvectors, Eigenvectors species, stress value and Scores from Ordination for envfit function
  Eigenvectors <- data.frame(vegan::scores(NMDS, choices=1:dimensions))
  EigenvectorsSpecies <- data.frame(NMDS$species)
  StressValue <- NMDS$stress
  
  # Save results
  results <- list(Eigenvectors=Eigenvectors, EigenvectorsSpecies=EigenvectorsSpecies, StressValue=StressValue, Ordination=NMDS)
  
  return(results)
}
