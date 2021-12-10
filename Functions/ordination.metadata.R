ordination.metadata <- function(metadata, ordination, permut = 1000, dimensions = c(1,2), p.adj.method = "BH", metadata.select){

  #' @title Fitting environmental data onto an ordination
  #'
  #' @description This function forms a wrapper around the function 'envfit' implemented in vegan. It fits environmental vectors or factors onto an ordination. The `significance' of fitted vectors or factors is assessed using permutations of environmental variables.
  #' Resulting p-values are corrected by the chosen adjustment method. Environmental vectors are scaled by their correlation values, consequentially the length of the resulting arrows is a measure for the correlation strength.
  #' The function returns a list with two elements, first a data table with fitted environmental vectors, second a list containing a data table with environmental factors and a vector with corresponding p-values.
  #'
  #' @section Required packages:
  #' 'vegan'
  #'
  #' @param metadata (\emph{data.frame/matrix}). A data.frame or matrix providing metadata. Column names must be factors, row names must be sample-IDs.
  #' @param ordination (\emph{ordination object}). An ordination object or other structure from which the ordination scores can be extracted (including a data frame or matrix of scores).
  #' @param permut (\emph{integer}). Number of permutations to check the significance of the fitted environmental variables. By default, permut = 1000.
  #' @param dimensions (\emph{integer}). A vector specifying two principal components used for fitting environmental data. By default, dimensions = c(1,2).
  #' @param p.adj.method (\emph{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}). Method, supported by 'p.adjust' implemented in stats, to correct p-values of the fitted metadata for multiple testing. By default, p.adj.method = "BH".
  #' @param metadata.select (\emph{integer/character}). Optional integer or characters.  If the argument is not defined all correlations are returned, else either a defined number of metadata vectors/factors (one \strong{integer} defining respectively a number of metadata factors and vectors) or only chosen factors(vectors) (a vector containing factor/vector names as \strong{characters}) are returned.
  
  #############################################################################################
  #################################### Check prerequisites ####################################

  #--------------------------- Check if needed arguments are valid ---------------------------#
  # Check if 'metadata' are provided
  if(missing(metadata) | !(class(metadata) %in% c("data.frame", "matrix"))){
    stop("Please provide a 'data.frame/matrix' with metadata.")
  }

  # Convert metadata to data.frame
  if(class(metadata) == "matrix"){ # class data = matrix
    metadata <- data.frame(metadata)
  } 

  # Check if 'ordination' are provided
  if(missing(ordination)){
    stop("Please provide an ordination object or other structure from which the ordination scores can be extracted.")
  }

  # Check if the argument 'permut' is valid
  if(!is.numeric(permut) | length(permut) != 1){
    stop("Please provide an integer defining the number of permutations for checking the significance of the fitted environmental variables.")
  } else if(permut%%1 != 0){
    stop("Please provide an integer defining the number of permutations for checking the significance of the fitted environmental variables.")
  }

  # Check if the argument 'dimensions' is valid
  if(!is.numeric(dimensions) | length(dimensions) != 2){
    stop("Please provide a vector containing two integers defining two principal components used for fitting environmental data.")
  } else if(any(dimensions%%1 != 0)){
    stop("Please provide a vector containing two integers defining two principal components used for fitting environmental data.")
  } else if(any(dimensions > ncol(vegan::scores(ordination)))){
    stop("One or both of the chosen principal components could not be found in the ordination results.")
  }

  # Check if the argument 'p.adj.method' is valid
  if(!p.adj.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")){
    stop("Wrong value for the argument 'p.adj.method' is provided. Please choose between 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'.")
  }

  #-------------------------- Check if optional arguments are valid --------------------------#

  # Check if the argument 'metadata.select' is valid
  if(!missing(metadata.select)){
    if(!(is.numeric(metadata.select) | is.character(metadata.select))){
      stop("Please provide an integer or characters, defining either the number of metadata (integer) with the highest correlations/smallest p-values or specific metadata factors (characters). Optional, if the argument is not defined, all correlations are returned.")
    } else if(is.numeric(metadata.select)){
      if(length(metadata.select) != 1 | metadata.select%%1 != 0){
        stop("Please provide an integer or characters, defining either the number of metadata (integer) with the highest correlations/ samllest p-values or specific metadata factors (characters). Optional, if the argument is not defined, all correlations are returned.")
      }
    }else if(!all(metadata.select %in% colnames(metadata))){
      stop("Please check the argument 'metadata.select'. All defined metadata factors must be listed in the column names of the metadata.")
    }
  } else {
    metadata.select <- NULL
  }

  #############################################################################################
  ################################ Fit metadata onto ordination ###############################

  # Fit environmental vectors or factors onto an ordination.
  ef <- vegan::envfit(vegan::scores(ordination, choices=dimensions), metadata, perm=permut)

  # Correct for multiple testing
  pvals.vectors <- p.adjust(ef$vectors$pvals, method=p.adj.method)
  pvals.factors <- p.adjust(ef$factors$pvals, method=p.adj.method)

  #--------------------------- Handling numeric metadata (vectors) ---------------------------#

  # Sort p-values in ascending order of the correlation values
  pvals.vectors=pvals.vectors[match(names(sort(ef$vectors$r, decreasing = T)), names(pvals.vectors))]

  # Only keep the chosen amount of top metdata or the selected metadata or all data if argument is not defined
  if(is.numeric(metadata.select)){
    if(length(pvals.vectors) > metadata.select){
      pvals.vectors <- pvals.vectors[1:metadata.select]
    }
  } else if(is.character(metadata.select)){
    pvals.vectors <- pvals.vectors[names(pvals.vectors) %in% metadata.select]
  }


  if(length(pvals.vectors) > 0){

    ## Extract arrows for numeric metadata
    metadata.arrows=data.frame(ef$vectors$arrows[match(names(pvals.vectors), rownames(ef$vectors$arrows)), ])

    ## If just one metdata is chose, for some reasons the matrix 'ef.arrows' has an different format
    ## Therefore one adjust its format as it should be (metadata as rows, dimensions as column)
    if(length(pvals.vectors) == 1){
      metadata.arrows=data.frame(t(metadata.arrows))
      row.names(metadata.arrows) <- names(pvals.vectors)
    }

    ## Change colnames of table
    colnames(metadata.arrows) <- c("x1", "y1")

    ## In plot the arrow values are scaled by their correlation (square root of the column r2) so that
    ## "weak" predictors have shorter arrows than "strong" predictors.
    ## For that, extract correlation values
    correlations=ef$vectors$r[match(names(pvals.vectors), names(ef$vectors$r))]

    ## It is possible that metadata.select is larger than the metadata in the file
    ## Correct to maximum possible metadata.select
    if(is.numeric(metadata.select)){
      if(length(metadata.arrows[,1]) < metadata.select){
        metadata.select <- length(metadata.arrows[,1])
      }
    }

    ## Scale values with corrlation value and, if wanted, with a chosen value (metadata.arrow.factor)
    metadata.arrows$x1 <- metadata.arrows$x1*correlations
    metadata.arrows$y1 <- metadata.arrows$y1*correlations

    ## Set labels and x0, y0
    metadata.arrows$Labels <- rownames(metadata.arrows)
    metadata.arrows$x0 <- 0
    metadata.arrows$y0 <- 0

    ## Assign p-values
    metadata.arrows$p.values <- pvals.vectors
  } else { # no vector found
    metadata.arrows <- NULL
  }

  #-------------------------- Handling categoric metadata (factors) --------------------------#

  # Sort factors by p-values
  pvals.factors=sort(pvals.factors)

  # Only keep the chosen amount of top metdata or the selected metadata or all data if argument is not defined
  if(is.numeric(metadata.select)){
    if(length(pvals.factors) > metadata.select){
      pvals.factors <- pvals.factors[1:metadata.select]
    }
  } else if(is.character(metadata.select)){
    pvals.factors <- pvals.factors[names(pvals.factors) %in% metadata.select]
  }

  if(length(pvals.factors) > 0){

    ## Extract factor data
    centroids <- data.frame(ef$factors$centroids[grep(paste(names(pvals.factors), collapse = "|"), rownames(ef$factors$centroids)), ])
    colnames(centroids) <- c("x", "y")

    ## Create a column for factors and names of factors respectively
    centroids$Factors <- rownames(centroids)
    centroids$Labels <- rownames(centroids)
    for(factor in names(pvals.factors)){
      centroids$Factors <- gsub(paste0(factor, ".*"), factor, centroids$Factors)
      centroids$Labels <- gsub(factor, "", centroids$Labels)
    }

    centroids <- list(centroids=centroids, `p-values`=pvals.factors)
  } else {
    centroids <- NULL
  }

  return(list(`Environmental vectors`=metadata.arrows, `Environmental factors`=centroids))
}
