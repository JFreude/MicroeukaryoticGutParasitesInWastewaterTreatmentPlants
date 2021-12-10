line.area.plot.transformation <- function(counts, metadata, taxonomy, x.axis, x.axis.order, fill.factor, fill.levels, split.fill.factor, split.number, fill.split.order = FALSE){

  #' @title Formatting community data for area and line plots
  #'
  #' @description This function converts community data into a format that allows creating area or line plots with ggplot. For this purpose mean values and standard deviations of all taxa are calculated across samples.
  #' The function returns a data set with the specified fill, (split), x-axis, mean values, and standard deviation, all listed in one column respectively. Also, a new column, \strong{Fill}, has been added, which assigns a letter to all fill and split values. This line must be specified as fill values for the area plots.
  #'
  #' @section Packages the function bases on:
  #' 'dplyr', 'reshape2', 'tidyr', 'plyr'
  #'
  #' @param counts (\emph{data.frame, matrix}). A data.frame or matrix providing count data. Column names must be sample-IDs, row names taxa-IDs. 
  #' @param metadata (\emph{data.frame, matrix}). A data.frame or matrix providing metadata. Column names must be factors, row names sample-IDs.
  #' @param taxonomy (\emph{data.frame, matrix}). A data.frame or matrix providing taxonomy data. Column names must be taxonomic ranks, row names taxa-IDs.
  #' @param x.axis (\emph{character}). A character specifying one column name of the metadata, with either characters/factors or integers. The values of the selected columns are plotted on the x-axis of the resulting area or line plot.
  #' @param x.axis.order (\emph{character}). It's an optional vector. If the argument \strong{x.axis} specifies a column containing characters/factors, this argument can be used to determine the order in which the values are plotted on the x-axis. All values of the column specified with the argument \strong{x.axis} must be found in this vector.
  #' @param fill.factor (\emph{character}). A character specifying one column name of the taxonomy data. The values of the selected columns are the fill values of the resulting area or line plot.
  #' @param fill.levels (\emph{integer}). An optional integer specifying the number of fill values. The total number of reads per taxon is calculated. The specified number of taxa with the highest reads is preserved while the other taxa are summed to one binned taxon. If the argument is not defined no taxa are summed up.
  #' @param split.fill.factor (\emph{character}). An optional character that specifies one column name of the taxonomy data. When this argument is set, it defines new fill values depending on the values defined by the argument \strong{fill.factor}. For example, if this argument is set to "genus" and \strong{fill.factor} is set to "Trait",  the resulting data frame would provide all mean values for genera per traits.
  #' @param split.number (\emph{integer}). An optional integer specifying the number of taxa per fill values. The specified number of taxa with the highest reads is preserved while the other taxa are summed to one binned taxon. If the argument is not defined no taxa are summed up.
  #' @param fill.split.order (\emph{logical}). If \emph{TRUE}, fill values are arranged in a way that taxa with high abundance at the first and low abundance at the last x-axis value are placed at the bottom, while taxa with opposite curves are placed at the top of the area plot. If \emph{FALSE}, taxa are arranged alphabetically. By default, fill.split.order = \emph{FALSE}.
  
  #############################################################################################
  #################################### Check prerequisites ####################################

  #----------------- Check if all needed functions and packages are available ----------------#
  if(!"dplyr" %in% rownames(installed.packages())){
    stop("Please install the package 'dplyr'.")
  }
  
  if(!"reshape2" %in% rownames(installed.packages())){
    stop("Please install the package 'reshape2'.")
  }
  
  if(!"tidyr" %in% rownames(installed.packages())){
    stop("Please install the package 'tidyr'.")
  }
  
  if(!"plyr" %in% rownames(installed.packages())){
    stop("Please install the package 'plyr'.")
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
  
  # Check if 'taxonomy' data is provided 
  if(missing(taxonomy) | !(class(taxonomy) %in% c("data.frame", "matrix"))){
    stop("Please provide taxonomy data formatted as data.frame or matrix")
  } else if(class(taxonomy) == "matrix"){
    taxonomy <- as.data.frame(taxonomy)
  } 
  
  # Check if corresponding rows and columns match
  if(all.equal(rownames(counts), rownames(taxonomy)) != TRUE){
    stop("Please check abundance and taxonomy data, they row names must match.")
  }
  if(all.equal(colnames(counts), rownames(metadata)) != TRUE){
    stop("Please check abundance and metadata, the column names of the abundance data must match the row names of the metadata. ")
  }

  # Check if the value of the argument 'x.axis' is valid
  if(!x.axis %in% colnames(metadata)){
    stop("Wrong value for the argument 'x.axis' is provided. Please select one column name of the metadata.")
  } else if(!class(metadata[,x.axis]) %in% c("factor", "character", "integer")){
    stop("Wrong value for the argument 'x.axis' is provided. The chosen column of the metadata must contain integers, characters or factors.")
  }
  if(class(metadata[,x.axis]) == "factor"){
    metadata[,x.axis] <- as.character(metadata[,x.axis])
  }

  # Check if the value of the argument 'x.axis.order' is valid
  if(missing(x.axis.order)){
    x.axis.order <- NULL
  } else {
    if(!all(as.character(unique(metadata[,x.axis])) %in% x.axis.order)){
      stop("Wrong value for the argument 'x.axis.order' is provided. All levels given for the selected fill factor must be found in the provided vector.")
    }
  }

  # Check if the value of the argument 'fill.factor' is valid
  if(!fill.factor %in% colnames(taxonomy)){
    stop("Wrong value for the argument 'fill.factor' is provided. Please select one column names of the taxonomy data.")
  }

  #-------------------------- Check if optional arguments are valid --------------------------#

  # Check if the value of the argument 'fill.levels' is valid
  if(!missing(fill.levels)){
    if(!is.numeric(fill.levels) | length(fill.levels) != 1){
      stop("Wrong value for the argument 'numb.taxa' is provided. Please define an integer specifying how many levels of the fill factor shall be shown or leave the argument empty.")
    } else if(fill.levels%%1 != 0){
      stop("Wrong value for the argument 'numb.taxa' is provided. Please define an integer specifying how many levels of the fill factor shall be shown or leave the argument empty.")
    }
  } else {# argument numb.taxa is missing
    fill.levels <- NULL
  }

  # Check if the value of the argument 'split.fill.factor' is valid
  if(missing(split.fill.factor)){
    split.fill.factor <- NULL
  } else if(!split.fill.factor %in% colnames(taxonomy)){
    stop("Wrong value for the argument 'split.fill.factor' is provided. Please select one column names of the taxonomy data.")
  }

  # Check if the value of the argument 'split.number' is valid
  if(!missing(split.number)){
    if(!is.numeric(split.number) | length(split.number) != 1){
      stop("Wrong value for the argument 'split.number' is provided. Please define an integer specifying how many levels of the split factor shall be shown or leave the argument empty.")
    } else if(split.number%%1 != 0){
      stop("Wrong value for the argument 'split.number' is provided. Please define an integer specifying how many levels of the split factor shall be shown or leave the argument empty.")
    }
  } else {# argument numb.taxa is missing
    split.number <- NULL
  }

  # Check if the value of the argument 'fill.split.order' is valid
  if(!fill.split.order %in% c(TRUE, FALSE)){
    stop("Wrong value for the argument 'fill.split.order' is provided. Please provide a logical (TRUE/FALSE).")
  }
  
  #############################################################################################
  ################################### NO 'split.fill.factor'###################################

  if(is.null(split.fill.factor)){
    ## Subset taxonomy data to get only used columns
    tax <- dplyr::select(taxonomy, all_of(fill.factor))

    ## Rename columns
    colnames(tax) <- c("Fill_col")

    ## Convert elements in columns to characters
    tax$Fill_col <- as.character(tax$Fill_col)

    ## Combine tax and abundance data
    data <- cbind(tax, counts)

    ## Sum data
    data <- aggregate(. ~ Fill_col, data, sum)

    #--------------------------------------- fill factor ---------------------------------------#

    ## Calculate the sum of reads per taxa per fill factor and sort data by sums
    data$sum <- rowSums(data[,which(lapply(data, class) == "numeric")])

    ## Sort by sums, highst to lowest value, to be able to select levels with highest counts
    data  <- data[order(-data$sum),]

    if(is.null(fill.levels)){
      ### Keep all taxa
      data <- data[ ,!colnames(data) %in% "sum"]
    } else {

      ### Check if there are more than wanted fill levels
      if(nrow(data) > fill.levels){

        #### Keep the most abundant taxa
        dat_sub <- data[seq(fill.levels),]

        #### Sum others up
        Others <- colSums(data[seq(fill.levels+1, nrow(data)),which(lapply(data, class) == "numeric")])

        #### Merge dat_sub and Other
        #### The assignment of the sum starts with z to enable sorting as wanted, later the name is changed back
        Others <- data.frame(Fill_col = "ZZZZZZZZ_Others", t(Others))

        #### Merge both dataframes
        data <- rbind(dat_sub, Others)
        data <- data[ ,!colnames(data) %in% "sum"]

      } else {
        ### Keep all taxa
        data <- data[ ,!colnames(data) %in% "sum"]
      }
    }

    ## Reshap data, each phylogenetic group, each treatment and the counts should all be in separate columns
    data <- reshape2::melt(data,id=c("Fill_col"))

    ## Merge with metadata
    metadata$variable <- rownames(metadata)
    data <- merge(data, metadata, by="variable", all = T)

    ## Calculate the mean and sd of the replicates (locations)
    sd <- plyr::ddply(data, c("Fill_col", x.axis),dplyr::summarise,sd=sd(value))$sd
    data <- plyr::ddply(data, c("Fill_col", x.axis),dplyr::summarise, Mean=mean(value))
    data$`Standard Deviation` <- sd

    ## Create new column with integers for x.axis
    if(class(metadata[,x.axis]) != "integer"){
      if(is.null(x.axis.order)){
        data <- data[order(factor(data[,x.axis])), ]
        data$X <- sort(rep(seq(unique(metadata[,x.axis])), nrow(data)/length(unique(metadata[,x.axis]))))
      } else {
        data <- data[order(factor(data[,x.axis], levels = x.axis.order)), ]
        data$X <- sort(rep(seq(x.axis.order), nrow(data)/length(x.axis.order)))
      }
    }

    ## Sort the dataset x-, and fill values as wanted, create new column with ascending letters per fill
    data <- data[order(factor(data$Fill_col) , data[x.axis]), ]
    data$Fill <- sort(rep(letters[seq(length(data$Fill_col)/length(unique(metadata[,x.axis])))],
                          nrow(data)/length(unique(data$Fill_col))))

    #------------------------------------- fill split order ------------------------------------#
    ## Split levels in decreasing order along x.axis
    if(fill.split.order){

      ### Spread values for x-axis
      data_spread <- data[, c("X","Fill","Mean")]
      data_spread <- tidyr::spread(data_spread, X, Mean)

      ### Set column Fill as rownames and delet it
      rownames(data_spread) <- data_spread$Fill
      data_spread <- data_spread[,!colnames(data_spread) %in% "Fill"]

      ### Order columns
      #for(i in  rev(colnames(data_spread))){
      #  data_spread <- data_spread[order(-data_spread[,i]),]
      #}
      data_spread$Order <- data_spread[,1]-data_spread[,ncol(data_spread)]
      data_spread <- data_spread[order(data_spread$Order),]

      ### Order data_to_sort
      data <- data[order(factor(data$Fill, levels = rownames(data_spread))),]

      ### Get others in last position
      Others <- sort(unique(data$Fill), decreasing = T)[1]
      data <- rbind(data[data$Fill != Others,], data[data$Fill == Others,])

      ### Assign new fill values
      data$Fill <- sort(data$Fill)
    }
    #-------------------------------------------------------------------------------------------#

    ## Renaming "ZZZZZZZZ_"
    data$Fill_col <- gsub("ZZZZZZZZ_", "", data$Fill_col)

    ## Renaming Fill_col
    colnames(data)[colnames(data)=="Fill_col"] <- fill.factor

  #############################################################################################
  ##################################### 'split.fill.factor'####################################
  } else {
    ## Subset taxonomy data to get only used columns
    tax <-  dplyr::select(taxonomy, fill.factor, split.fill.factor)

    ## Rename columns
    colnames(tax) <- c("Fill_col", "Split_col")

    ## Convert elements in columns to charcters
    tax[, ] <- lapply(tax[, ], as.character)

    ## Combine tax and abundance data
    data <- cbind(tax, counts)

    ## Sum data
    data <- aggregate(. ~ Fill_col+Split_col, data, sum)
    #------------------------------------- split fill factor -----------------------------------#

    ## Create separate data frames for fill factors
    data <- lapply(dplyr::group_split(data, Fill_col), data.frame)

    ## Create open list
    dataset <- list()

    ## Iterate over data subsets
    for(level in seq(data)){
      dat <- data[[level]]

      ### Calculate the sum of reads per taxa per fill factor and sort data by sums
      dat$sum <- rowSums(dat[,which(lapply(dat, class) == "numeric")])

      ### Sort by sums, highst to lowest value, to be able to select levels with highest counts
      dat  <- dat[order(-dat$sum),]

      if(is.null(split.number)){
        #### Keep all taxa
        dat <- dat[ ,!colnames(dat) %in% "sum"]
      } else {

        #### Check if there are more than wanted split levels
        if(nrow(dat) > split.number){

          ##### Keep the most abundant taxa
          dat_sub <- dat[seq(split.number),]

          ##### Sum others up
          Others <- colSums(dat[seq(split.number+1, nrow(dat)),which(lapply(dat, class) == "numeric")])

          ##### Merge dat_sub and Other
          ##### The assignment of the sum starts with z to enable sorting as wanted, later the name is changed back
          Others <- data.frame(Fill_col = dat$Fill_col[1], Split_col = paste("ZZZZZZZZ_Other", dat$Fill_col[1]), t(Others))

          ##### Merge both dataframes
          dat <- rbind(dat_sub, Others)
          dat <- dat[ ,!colnames(dat) %in% "sum"]

        } else {
          ##### Keep all taxa
          dat <- dat[ ,!colnames(dat) %in% "sum"]
        }
      }

      dataset[[level]] <- dat
    }# for loop

    ## Merge all data
    data <- do.call(rbind, dataset)

    #---------------------------------------- fill level ---------------------------------------#

    if(!is.null(fill.levels)){

      ### Sum values up by fill
      data_fill <- data[,!colnames(data) %in% "Split_col"]
      data_fill <- aggregate(. ~ Fill_col, data_fill, sum)

      ### Calculate the sum per fill level
      data_fill$sum <- rowSums(data_fill[ ,which(lapply(data_fill, class) == "numeric")])

      ### Sort by sum, highst to lowest value, to be able to select levels with highest reads
      data_fill <- data_fill[order(-data_fill$sum), ]

      ### Keep the 15 most abundant taxa
      if(nrow(data_fill) > fill.levels){

        #### Get names of wanted levels
        fill_names <- data_fill[seq(1,fill.levels),"Fill_col"]

        #### Keep the wanted levels
        data_fill <- data[data[,"Fill_col"] %in% fill_names, ]

        #### Sum up other taxa
        #### The assignment of the sum starts with z to enable sorting as wanted, later the name is changed back
        Others <- data[!data[,"Fill_col"] %in% fill_names, ]
        Others$Fill_col <- "ZZZZZZZZ_Others"
        Others <- aggregate(. ~ Fill_col+Split_col, Others, sum)

        #### Merge both dataframes
        data <- rbind(data_fill, Others)
      }
    }

    ## Reshap data, each phylogenetic group, each treatment and the abundaces should all be in seperate columns
    data <- reshape2::melt(data,id=c("Fill_col","Split_col"))

    ## Merge with metadata
    metadata$variable <- rownames(metadata)
    data <- merge(data, metadata, by="variable", all = T)

    ## Calculate the mean and sd of the replicates (locations)
    sd <- plyr::ddply(data, c("Fill_col", "Split_col", x.axis),dplyr::summarise,sd=sd(value))$sd
    data <- plyr::ddply(data, c("Fill_col", "Split_col", x.axis),dplyr::summarise, Mean=mean(value))
    data$`Standard Deviation` <- sd

    ## Create new column with integers for x.axis
    if(class(metadata[,x.axis]) != "integer"){
      if(is.null(x.axis.order)){
        data <- data[order(factor(data[,x.axis])), ]
        data$X <- sort(rep(seq(unique(metadata[,x.axis])), nrow(data)/length(unique(metadata[,x.axis]))))
      } else {
        data <- data[order(factor(data[,x.axis], levels = x.axis.order)), ]
        data$X <- sort(rep(seq(x.axis.order), nrow(data)/length(x.axis.order)))
      }
    }

    ## Sort the dataset x-, and fill values as wanted, create new column with ascending letters per fill
    data <- data[order(factor(data$Fill_col), data$Split_col , data[x.axis]), ]
    data$Fill <- sort(rep(letters[seq(length(data$Split_col)/length(unique(metadata[,x.axis])))],
      nrow(data)/length(unique(data$Split_col))))

    #------------------------------------- fill split order ------------------------------------#
    ## Split levels in decreasing order along x.axis
    if(fill.split.order){
      ### Create separate data frames for fill factors
      data_sort <- lapply(dplyr::group_split(data, Fill_col), data.frame)

      ### Open list
      data_sorted <- list()

      ### Iterate over data subsets
      for(level in seq(data_sort)){
        data_to_sort <- data_sort[[level]]

        #### Change colnames for sd
        colnames(data_to_sort)[colnames(data_to_sort) == "Standard.Deviation"] <- "Standard Deviation"

        #### Spread values for x-axis
        data_spread <- data_to_sort[, c("X","Fill","Mean")]
        data_spread <- tidyr::spread(data_spread, X, Mean)

        #### Set column Fill as rownames and delet it
        rownames(data_spread) <- data_spread$Fill
        data_spread <- data_spread[,!colnames(data_spread) %in% "Fill"]

        #### Order columns
        #for(i in  rev(colnames(data_spread))){
        #  data_spread <- data_spread[order(-data_spread[,i]),]
        #}
        data_spread$Order <- data_spread[,1]-data_spread[,ncol(data_spread)]
        data_spread <- data_spread[order(data_spread$Order),]

        #### Order data_to_sort
        data_to_sort <- data_to_sort[order(factor(data_to_sort$Fill, levels = rownames(data_spread))),]

        #### Get others in last position
        Others <- sort(unique(data_to_sort$Fill), decreasing = T)[1]
        data_to_sort <- rbind(data_to_sort[data_to_sort$Fill != Others,], data_to_sort[data_to_sort$Fill == Others,])

        #### Assign new fill values
        data_to_sort$Fill <- sort(data_to_sort$Fill)

        #### Save results in list
        data_sorted[[level]] <- data_to_sort
      }

      ### Merge lists
      data <- do.call(rbind, data_sorted)
    }
    #-------------------------------------------------------------------------------------------#

    ## Renaming "ZZZZZZZZ_"
    data$Split_col <- gsub("ZZZZZZZZ_", "", data$Split_col)
    data$Fill_col <- gsub("ZZZZZZZZ_", "", data$Fill_col)

    ## Renaming Split_col, Fill_col
    colnames(data)[colnames(data)=="Split_col"] <- split.fill.factor
    colnames(data)[colnames(data)=="Fill_col"] <- fill.factor

  }

  return(data)
}
