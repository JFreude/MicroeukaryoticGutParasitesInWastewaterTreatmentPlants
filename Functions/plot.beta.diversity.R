plot.beta.diversity <- function(plot.data, wilcoxon.pvalues=TRUE, paired=FALSE, colored=TRUE, shaped=FALSE){  
  
  #' @title Plot beta diversity across groups
  #' @description The function allows to compare beta diversities across groups by plotting them as box-plot. 
  #' Optionally, a two-sided Wilcoxon test is conducted using the function stat_compare_means implemented in ggpubr, the significant results (p < 0.05) are shown in the plot.
  #' 
  #' @param plot.data (\emph{matrix/dataframe}). A matrix/dataframe with beta diverity indices per group as columns and samples as row. Additionally, if desired, columns with group assignments for colours and shapes can also be included in the data set. According to these groups the samples are colored and shaped.
  #' @param wilcoxon.pvalues (\emph{logical}). If \emph{TRUE}, significant Wilcoxon p-values are displayed on the box plot using function stat_compare_means in R package ggpubr.
  #' @param paired (\emph{logical}). If \emph{TRUE} a paired, if \emph{FALSE} an unpaired wilcoxon test is conducted. By default an unpaired wilcoxon test is conducted. 
  #' @param colored (\emph{logical/character}). If \emph{TRUE}, the samples are colored by the samples, if \emph{FALSE}, they are not colored. Alternatively, a \emph{character}, one column name of the \strong{plot.data}, can be specified. By default the samples are colored by the samples. 
  #' @param shaped (\emph{logical/character}). If \emph{TRUE}, the samples are shaped by the samples, if \emph{FALSE} they are not shaped. Alternatively, a \emph{character}, one column name of the \strong{plot.data}, can be specified. By default the samples are not shaped. 
  #' 
  #' @section Packages the function bases on
  #' 'ggplot2', 'ggpubr', 'reshape2'
  #' 
  
  #############################################################################################
  #################################### Check prerequisites ####################################
  
  #----------------- Check if all needed functions and packages are available ----------------#
  
  if(!"ggplot2" %in% rownames(installed.packages())){
    stop("Please install the package 'ggplot2'.")
  } 
  
  if(!"ggpubr" %in% rownames(installed.packages())){
    stop("Please install the package 'ggpubr'.")
  }
  
  if(!"reshape2" %in% rownames(installed.packages())){
    stop("Please install the package 'reshape2'.")
  } 
  
  #--------------------------- Check if needed arguments are valid ---------------------------#
  
  # Check if 'plot.data' are provided 
  if(missing(plot.data) | !(class(plot.data) %in% c("data.frame", "matrix"))){
    stop("Please provide a 'data.frame' or 'matrix' with diverity indices per group as columns and samples as row.")
  } else if(class(plot.data) == "matrix") {
    plot.data <- as.data.frame(plot.data)
  }
  
  # Check if the value of the argument 'wilcoxon.pvalues' is valid
  if(!wilcoxon.pvalues %in% c(TRUE, FALSE)){
    stop("Wrong value for the argument 'wilcoxon.pvalues' is provided. Please provide a logical (TRUE/FALSE).")
  } 
  
  # Check if the value of the argument 'paired' is valid
  if(!paired %in% c(TRUE, FALSE)){
    stop("Wrong value for the argument 'paired' is provided. Please provide a logical (TRUE/FALSE).")
  } 
  
  
  # Check if the value of the argument 'colored' is valid
  if(!colored %in% c(TRUE, FALSE, colnames(plot.data))){
    stop("Wrong value for the argument 'colored' is provided. Please provide a logical (TRUE/FALSE) or a coloumn name of plot.data.")
  } 
  
  # Check if the value of the argument 'colored' is valid
  if(!shaped %in% c(TRUE, FALSE, colnames(plot.data))){
    stop("Wrong value for the argument 'shaped' is provided. Please provide a logical (TRUE/FALSE) or a coloumn name of plot.data.")
  } 

  #############################################################################################
  ####################################### Diversity Plot ######################################
  
  #--------------------------- Calculate p-values by Wilcoxon test  --------------------------#
  
  # Display box plot with p-values
  if(wilcoxon.pvalues){
    
    # Subset data, keep only diversity indices
    wilcoxon_data <- plot.data[,!colnames(plot.data) %in% c(colored, shaped)]
    
    # Get significant paires
    combinations=list()
    units=colnames(wilcoxon_data)
    
    # Iterates each column of 'wilcoxon_data'
    for(index1 in 1:(length(units)-1)){
      for(index2 in (index1+1):length(units)){
        unit1=units[index1]
        unit2=units[index2]
        
        # Two-sided, paired Wilcoxon test
        if(paired){
          w.out=wilcox.test(wilcoxon_data[,index1],wilcoxon_data[,index2], paired = TRUE)
          # Two-sided, unpaired Wilcoxon test
        }else{
          w.out=wilcox.test(wilcoxon_data[,index1],wilcoxon_data[,index2], paired = FALSE)
          print(w.out)
          print(wilcox.test(wilcoxon_data[!is.na(wilcoxon_data[,index1]),index1], wilcoxon_data[!is.na(wilcoxon_data[,index2]),index2], paired = FALSE))
        }
        
        if(w.out$p.value<0.05){
          combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2)
        }
      }
    }
    
    print(paste("Number of significant differences across groups:",length(combinations)))
  }
  
  #-------------------------------------- Plot diversity -------------------------------------#  
  
  # Transform data.frame for plotting
  df_melt <- suppressMessages(reshape2::melt(as.matrix(plot.data[,!colnames(plot.data) %in% c(colored, shaped)])))
  colnames(df_melt) <- c("Samples", "Groups", "value")
  
  # Remove NAs
  df_melt <- df_melt[!is.na(df_melt$value),] 
  
  # Add columns for colors and shapes
  if(!all(colored %in% c(TRUE, FALSE))){
    df_melt <- setNames(cbind(df_melt, as.character(plot.data[match(df_melt$Samples, rownames(plot.data)), colored])), 
                        c(colnames(df_melt), colored))
  }
  if(!all(shaped %in% c(TRUE, FALSE))){
    df_melt <- setNames(cbind(df_melt, as.character(plot.data[match(df_melt$Samples, rownames(plot.data)), shaped])), 
                        c(colnames(df_melt), shaped))
  }
  
  # Separate Sample Pairs
  # Check if Sample names have just one point in Sample name, otherwise splitting will not work
  if(length(unlist(strsplit(as.character(df_melt$Samples), "\\."))) == 2*nrow(df_melt)){
    #df_melt <- data.frame(df_melt %>% separate("Samples", c("Sample1", "Sample2"), "\\."))
    df_melt <- cbind(Sample1 = gsub("\\..*$", "", df_melt$Samples), Sample2 = gsub("^.*\\.", "", df_melt$Samples),
                     df_melt[,colnames(df_melt) != "Samples"])
  } else {
    stop("Sample pairs can not be seperated by dots, please make shur that only the samples are separeted by dot (i.e. Sample1.Sample2).")
  }

  
  # Basic box-plot  
  p <- ggplot2::ggplot(df_melt,  ggplot2::aes(x=Groups, y=value))+ ggplot2::geom_boxplot()
  
  # If requested, add p-values
  if(wilcoxon.pvalues){
    # Two-sided, paired Wilcoxon test
    if(paired){
      p <- p + ggpubr::stat_compare_means(comparisons=combinations, paired = TRUE)
      # Two-sided, unpaired Wilcoxon test
    }else{
      p <- p + ggpubr::stat_compare_means(comparisons=combinations, paired = FALSE)
    }
  }
  
  # Color and shape of Samples
  
  #------ colored and shaped are same classes
  if(class(colored) == class(shaped)){
    ##------ colored and shaped are logicals
    if(class(colored) == "logical" & colored == T){
      p <- p + ggplot2::geom_jitter(ggplot2::aes(color=Sample1, shape=Sample2),
                                    position=ggplot2::position_jitter(0.2)) +
        ggplot2::scale_shape_manual(values=seq(unique(df_melt$Sample2)))
      ##------ colored and shaped are characters
    } else {
      p <- p + ggplot2::geom_jitter(ggplot2::aes(color=df_melt[,colored], shape=df_melt[,shaped]),
                                    position=ggplot2::position_jitter(0.2)) + 
        ggplot2::scale_shape_manual(values=seq(unique(df_melt[,shaped]))) +
        ggplot2::labs(color = colored, shape = shaped)
    }
    #------ colored and shaped are different classes
  } else {
    ##------ colored is logical, shape is character
    if(class(colored) == "logical"){
      if(colored == T){
        p <- p + ggplot2::geom_jitter(ggplot2::aes(color=Sample1, shape=df_melt[,shaped]),
                                      position=ggplot2::position_jitter(0.2)) + 
          ggplot2::scale_shape_manual(values=seq(unique(df_melt[,shaped]))) +
          ggplot2::labs(shape = shaped)
      } else {
        p <- p + ggplot2::geom_jitter(ggplot2::aes(shape=df_melt[,shaped]),
                                      position=ggplot2::position_jitter(0.2)) + 
          ggplot2::scale_shape_manual(values=seq(unique(df_melt[,shaped]))) +
          ggplot2::labs(shape = shaped)
      }
      ##------ colored is character, shape is logical
    } else {
      if(shaped == T){
        p <- p + ggplot2::geom_jitter(ggplot2::aes(color=df_melt[,colored], shape=Sample2),
                                      position=ggplot2::position_jitter(0.2)) + 
          ggplot2::scale_shape_manual(values=seq(unique(df_melt$Sample2))) +
          ggplot2::labs(color = colored)
      } else {
        p <- p + ggplot2::geom_jitter(ggplot2::aes(color=df_melt[,colored]),
                                      position=ggplot2::position_jitter(0.2)) + 
          ggplot2::labs(color = colored)
      }
    }
  }
  
  return(p)
}
