calculate.beta.diversity <-function(counts, metadata, groups, subsample=TRUE, normalize=TRUE){
  
  #' @title Calculate Bray-Curtis dissimilarity across groups
  #' @description The function calculates the Bray-Curtis dissimilarity, using the function vegdist implemented in vegan, across groups. 
  #' It returns a data set with dissimilarity indices per group as columns and samples as row. 
  #' Groups can be subsampled to the same sample number. This is a random process following the output of the same run can differ.
  #' 
  #' @param counts (\emph{data.frame, matrix}). A data.frame or matrix providing count data. Column names must be sample-IDs, row names taxa-IDs. 
  #' @param metadata (\emph{data.frame, matrix}). A data.frame or matrix providing metadata. Column names must be factors, row names must be sample-IDs.
  #' @param groups (\emph{character}). A character specifying one column name of the metadata, with either characters/factors or integers. The values of the selected columns are plotted on the x-axis of the resulting box plots.
  #' @param subsample (\emph{logical}). If \emph{TRUE}, groups are randomly subsampled to have the same sample number in each group, else data are not subsampled. By default, data are subsampled. 
  #' @param normalize (\emph{logical}). If \emph{TRUE}, count data are normalized row-wise (relative countss are calculated). By default, counts are normalized.
  #' 
  #' @section Packages the function bases on
  #' 'vegan', 'reshape2'

  
  #############################################################################################
  #################################### Check prerequisites ####################################
  
  #----------------- Check if all needed functions and packages are available ----------------#
  if(!"vegan" %in% rownames(installed.packages())){
    stop("Please install the package 'vegan'.")
  }
  
  if(!"reshape2" %in% rownames(installed.packages())){
    stop("Please install the package 'reshape2'.")
  }
  
  #--------------------------- Check if needed arguments are valid ---------------------------#
  
  # Check if 'counts' are provided
  if(missing(counts) | !(class(counts) %in% c("data.frame", "matrix"))){
    stop("Please provide count data formatted as data.frame or matrix")
  } else if(!all(unlist(lapply(counts, is.numeric)))){
    stop("Please check your count data, some of the entries are not numeric")
  } else if(class(counts) == "matrix"){
    counts <- as.data.frame(counts)
  }
  
  # Check if 'metadata' are provided 
  if(missing(metadata) | !(class(metadata) %in% c("data.frame", "matrix"))){
    stop("Please provide metadata formatted as data.frame or matrix")
  } else if(class(metadata) == "matrix"){
    metadata <- as.data.frame(metadata)
  } 

  # Check if corresponding rows and columns match
  if(all.equal(colnames(counts), rownames(metadata)) != TRUE){
    stop("Please check count- and metadata, the column names of the count data must match the row names of the metadata.")
  }
  
  # Check if the value of the argument 'groups' is valid
  if(!groups %in% colnames(metadata)){
    stop("Wrong value for the argument 'groups' is provided. Please select one column name of the metadata.")
  } else if(!class(metadata[,groups]) %in% c("factor", "character", "integer")){
    stop("Wrong value for the argument 'groups' is provided. The chosen column of the metadata must contain integers, characters or factors.")
  }
  if(class(metadata[,groups]) == "factor"){
    metadata[,groups] <- as.character(metadata[,groups])
  }
  
  # Check if the value of the argument 'subsample' is valid
  if(!subsample %in% c(TRUE, FALSE)){
    stop("Wrong value for the argument 'subsample' is provided. Please provide a logical (TRUE/FALSE).")
  } else if(subsample == TRUE){
    if(length(unique(table(metadata[,groups]))) == 1){
      subsample <- FALSE
    } else {
      print(paste("Constraining sample number randomly to the same minimal group sample number of", min(table(metadata[,groups]))))
    }
  }
  
  # Check if the value of the argument 'normalize' is valid
  if(!normalize %in% c(TRUE, FALSE)){
    stop("Wrong value for the argument 'normalize' is provided. Please provide a logical (TRUE/FALSE).")
  } 
  if(normalize){
    # Normalize column-wise = calculate relative counts
    counts=sweep(x = counts, MARGIN = 2, STATS = colSums(counts), FUN = '/')
    print("counts were normalizede column-wise")
  } else {
    warning("Please be aware that the abundance data are not normalized.")
  }
  
  #############################################################################################
  ######################################## Bray-curtis ########################################
  
  #------------------------ Prepare data to calculate the Bray-curtis ------------------------#
  
  # Unique group vector
  unique.groups <- unique(metadata[,groups])
  
  # If group vector is numeric, sort it
  if(is.numeric(unique.groups)){
    unique.groups <- sort(unique.groups)
  }
  
  #------------------- Calculate Bray-curtis per each member of each group -------------------#
  
  # Open data frame for Bray-curtis indices per group and sample
  Sample=expand.grid(unique(rownames(metadata)),unique(rownames(metadata)))
  beta_per_group <- data.frame(Sample=paste(Sample$Var1, Sample$Var2, sep="."))
  
  # Iterates over each group of the group vector
  for(group in unique.groups){
    ## Open list for Bray-curtis indices per each member of one group
    beta_per_group_member <- c()
    
    ## Extract indices of all group members
    group.indices=which(metadata[,groups]==group)
    
    ## If specified, subsample groups randomly to same sample number per group
    if(subsample==TRUE){
      group.indices=sample(group.indices)[1:min(table(metadata[,groups]))]
    }
    
    ## Subsample count data to only group members
    group.data <- counts[,group.indices]
    group.data <- group.data[,order(colnames(group.data))]
    
    # Calculate Bray curtis per group member
    if(is.null(group.data)==FALSE && ncol(group.data) > 1){
      dist=as.matrix(vegan::vegdist(t(group.data),method="bray"))
      dist[lower.tri(dist, diag = TRUE)] <- NA 
      dist <- subset(reshape2::melt(dist, na.rm = TRUE))
      values <- as.vector(dist$value)
      names(values) <- paste(dist$Var1, dist$Var2, sep="_____")
    } else{
      warning(paste("Group",group,"has less than 2 samples: cannot compute beta diversity."))
    }
    
    ## Convert vector to data frame
    beta_per_group_member <- setNames(data.frame(values), as.character(group))
    
    ## Create Column Sample 1 and Sample 2
    beta_per_group_member$Sample1 <- gsub("_____.*", "", rownames(beta_per_group_member))
    beta_per_group_member$Sample2 <- gsub(".*_____", "", rownames(beta_per_group_member))
    
    ## Extract sample names of metadata, merge them to one column
    beta_per_group_member$Sample <- paste(beta_per_group_member$Sample1, beta_per_group_member$Sample2, sep = ".")
    beta_per_group_member$Sample1 <- NULL
    beta_per_group_member$Sample2 <- NULL
    
    ## Merge all groups
    beta_per_group <- merge(beta_per_group, beta_per_group_member, by="Sample", all = T)
    
  } # end group loop
  
  # Set sample names as rownames
  rownames(beta_per_group) <- beta_per_group$Sample
  beta_per_group$Sample <- NULL
  
  # Remove rows with only NAs
  beta_per_group <- beta_per_group[rowSums(!is.na(beta_per_group))!=0,]
  
  return(beta_per_group)
}