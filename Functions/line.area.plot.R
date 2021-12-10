line.area.plot <- function(plot.data, plot.type="area.smooth", x.axis, y.axis, fill.order, sd, negative.sd=TRUE, legend.names, color){
  
  #' @title Area and line plots with ggplot
  #'
  #' @description This function can visualize area or line plots with sharp or smooth edges. The line plots show the standard deviations of the respective mean values. For visualization 'ggplot' is used.
  #'
  #' @section ackages the function bases on:
  #' 'dplyr', 'ggplot2', 'RColorBrewer'
  #'
  #' @param plot.data (\emph{data.frame/matrix}). A data.frame/matrix providing x-axis values (x.axis), mean values of taxa per sample (y.axis), the corresponding standard deviations (sd) and fill values (fill.order), each per column. Additionally, legend names and colors can be specified.
  #' @param plot.type (\emph{(\emph{"area.edged", "area.smooth", "line.edged", "line.smooth"}). A character specifying if an area- or line plot is plotted. The lines and areas of both can be plotted with sharp (\strong{'area.edged', 'line.edged'}) or smooth (\strong{'area.smooth', 'line.smooth'}) edges.
  #' @param x.axis (\emph{character}). A character specifying one column name of \strong{plot.data}, with numeric/integer x-axis values. The values of the selected columns are plotted on the x-axis of the resulting area or line plot.
  #' @param y.axis (\emph{character}). A character specifying one column name of \strong{plot.data}, with numeric y-axis values (mean value). The values of the selected columns are plotted on the y-axis of the resulting area or line plot.
  #' @param fill.order (\emph{character}). A character specifying one column name of \strong{plot.data} with fill values. The values of the selected columns are used as fill values of the resulting area or line plot. They are arranged alphabetically in the plot. Therefore, one may consider using letters for 'fill.order'. The corresponding names for the legend can be specified in the parameter \strong{plot.data}.
  #' @param sd (\emph{character}). A character specifying one column name of \strong{plot.data}, with numeric standard deviations. The argument must be defined if the argument \strong{plot.type} is set to \emph{'line.edged'} or \emph{'line.smooth'}. The values of the selected columns are used as error values for the line plot.
  #' @param negative.sd (\emph{logical}). The argument specifies whether deviations from the mean < 0 (mean - standard deviations < 0) are shown in the line plots. If \strong{TRUE} deviations from the mean < 0 are shown, if \strong{FALSE} the error bands of the line plots range only to 0.
  #' @param legend.names (\emph{character}). An optional character specifying one column name of \strong{plot.data} with new legend names. Only one legend name must be specified per \strong{fill.order} value. The values replace the fill values specified in the parameter \strong{fill.order} in the legend.
  #' @param color (\emph{character}). An optional character specifying one column name of \strong{plot.data} with colors.Only one colour must be specified per \strong{fill.order} value.
  

  #############################################################################################
  #################################### Check prerequisites ####################################

  #----------------- Check if all needed functions and packages are available ----------------#
  if(!"dplyr" %in% rownames(installed.packages())){
    stop("Please install the package 'dplyr'.")
  }
  
  if(!"ggplot2" %in% rownames(installed.packages())){
    stop("Please install the package 'ggplot2'.")
  }
  
  if(!"RColorBrewer" %in% rownames(installed.packages())){
    stop("Please install the package 'RColorBrewer'.")
  }
  
  #--------------------------- Check if needed arguments are valid ---------------------------#

  # Check if 'plot.data' are provided
  if(missing(plot.data) | !(class(plot.data) %in% c("data.frame", "matrix"))){
    stop("Please provide a data.frame/matrix with x-axis values (x.axis), mean values of taxa per sample (y.axis), the corresponding standard deviations (sd) and fill values (fill.order).")
  } else if(class(plot.data) == "matrix"){
    plot.data <- data.frame(plot.data)
  }

  # Check if the value of the argument 'plot.type' is valid
  if(!plot.type %in% c("area.edged", "area.smooth", "line.edged", "line.smooth")){
    stop("Wrong value for the argument 'plot.type' is provided. Please choose between 'area.edged', 'line.edged', 'area.smooth' and 'line.smooth'.")
  }

  # Check if the value of the argument 'x.axis' is valid
  if(!x.axis %in% colnames(plot.data)){
    stop("Wrong value for the argument 'x.axis' is provided. Please select one column name of the plot data.")
  } else if(!class(plot.data[,x.axis]) %in% c("numeric", "integer")){
    stop("Wrong value for the argument 'x.axis' is provided. The chosen column of the plot data must contain integers or numerics.")
  }

  # Check if the value of the argument 'y.axis' is valid
  if(!y.axis %in% colnames(plot.data)){
    stop("Wrong value for the argument 'y.axis' is provided. Please select one column name of the plot data.")
  } else if(!class(plot.data[,y.axis]) == "numeric"){
    stop("Wrong value for the argument 'y.axis' is provided. The chosen column of the plot data must contain numerics.")
  }

  # Check if the value of the argument 'fill.order' is valid
  if(!fill.order %in% colnames(plot.data)){
    stop("Wrong value for the argument 'fill.order' is provided. Please select one column name of the plot data.")
  }

  # Check if the value of the argument 'sd' is valid
  if(plot.type %in% c("line.edged", "line.smooth")){
    if(!sd %in% colnames(plot.data)){
      stop("Wrong value for the argument 'sd' is provided. Please select one column name of the plot data.")
    } else if(!class(plot.data[,sd]) == "numeric"){
      stop("Wrong value for the argument 'sd' is provided. The chosen column of the plot data must contain numerics.")
    }
  }

  # Check if the value of the argument 'negative.sd' is valid
  if(class(negative.sd) != "logical"){
    stop("Wrong value for the argument 'negative.sd' is provided. Please define a logical.")
  } else if(!plot.type %in% c("line.edged", "line.smooth")){
    negative.sd = TRUE
  }

  #-------------------------- Check if optional arguments are valid --------------------------#

  # Check if the value of the argument 'legend.names' is valid
  if(!missing(legend.names)){
    if(!legend.names %in% colnames(plot.data)){
      stop("Wrong value for the argument 'legend.names' is provided. Please select one column name of the plot data.")
    }
  }

  # Check if the value of the argument 'color' is valid
  if(!missing(color)){
    if(!color %in% colnames(plot.data)){
      stop("Wrong value for the argument 'color' is provided. Please select one column name of the plot data.")
    }
  }

  #############################################################################################
  ########################################### Plots ###########################################

  #--------------------------------- Data frame for plotting ---------------------------------#

  if(plot.type %in% c("line.edged", "line.smooth")){
    plot_data <- plot.data[,c(x.axis, y.axis, fill.order, sd)]
    colnames(plot_data) <- c("x","y","fill.order","sd")
  } else {
    plot_data <- plot.data[,c(x.axis, y.axis, fill.order)]
    colnames(plot_data) <- c("x","y","fill.order")
  }

  # Add column for colors
  if(!missing(legend.names)){
    plot_data$Fill <- plot.data[,legend.names]
  } else {
    plot_data$Fill <- plot.data[,fill.order]
  }

  # Add column for colors
  if(!missing(color)){
    plot_data$Col <- plot.data[,color]
  } else {
    cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(plot_data[,"fill.order"])))
    plot_data$Col <- sort(rep(cols,nrow(plot_data)/length(unique(plot_data[,"fill.order"]))))
  }

  # Convert factors to character, sort by 'fill.order'
  plot_data <- dplyr::mutate_if(plot_data, is.factor, as.character)
  plot_data <- plot_data[order(plot_data$fill.order),]

  ############################################################ NEW
  # No negative sd
  plot_data$sd_min <- plot_data$sd
  if(plot.type %in% c("line.edged", "line.smooth") & !negative.sd){
    negative_values <- plot_data$y - plot_data$sd < 0
    plot_data[negative_values,"sd_min"] <- plot_data[negative_values,"y"]
  }

  #------------------------------------- Edged area plot -------------------------------------#

  if(plot.type == "area.edged"){
    plot <- ggplot2::ggplot(NULL, ggplot2::aes(x=x, y=y, fill=fill.order)) +
      ggplot2::geom_area(data = plot_data, position = "stack") +
      ggplot2::geom_line(data = plot_data[plot_data$fill.order != plot_data$fill.order[1],], position = "stack", size = 0.5)
  }

  #------------------------------------- Smooth area plot ------------------------------------#

  if(plot.type == "area.smooth") {
    plot <- ggplot2::ggplot(NULL, ggplot2::aes(x=x, y=y, fill=fill.order)) +
      ggplot2::stat_smooth(data = plot_data,
                  geom = 'area', method = 'loess', alpha = 1, position = "stack") +
      ggplot2::stat_smooth(data = plot_data[plot_data$fill.order != plot_data$fill.order[1],], method = "loess", position = "stack", se = FALSE,
                  col="black", size = 0.5, show.legend=FALSE)
  }

  #------------------------------------- Edged line plot -------------------------------------#

  if(plot.type == "line.edged"){
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x=x, y=y, color=fill.order, fill=fill.order)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = y - sd_min, ymax = y + sd))
  }

  #------------------------------------- Smooth line plot ------------------------------------#

  if(plot.type == "line.smooth") {

    ## Add ymin and ymax to plot_data
    plot_data$ymin = plot_data$y - plot_data$sd_min
    plot_data$ymax = plot_data$y + plot_data$sd

    ## Create smooth lines for y, ymin, ymax
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(fill=fill.order)) +
      ggplot2::stat_smooth(ggplot2::aes(x = x, y = y, colour = "mean"), method = 'loess', se = FALSE) +
      ggplot2::stat_smooth(ggplot2::aes(x = x, y = ymin, colour = "min"), method = 'loess', se = FALSE) +
      ggplot2::stat_smooth(ggplot2::aes(x = x, y = ymax, colour = "max"), method = 'loess', se = FALSE)

    p <- ggplot2::ggplot_build(p)

    ## Extract values from plot
    y = setNames( p$data[[1]][,c("x","y","group")], c("x", "y", "fill.order") )
    ymin = setNames( p$data[[2]][,c("x","y","group")], c("x", "ymin", "fill.order") )
    ymax = setNames( p$data[[3]][,c("x","y","group")], c("x", "ymax", "fill.order") )

    # Create new table
    plot_data_new <- merge(merge(y,ymin, by = c("x","fill.order"), all = TRUE), ymax, by = c("x","fill.order"), all = TRUE)
    plot_data_new$Fill <- unique(plot_data$Fill)[plot_data_new$fill.order]
    plot_data_new$fill.order <- letters[plot_data_new$fill.order]
    plot_data_new <- merge(plot_data_new, plot_data[,c("Col", "Fill")], by="Fill", all = T)

    # Convert factors to character, sort by 'fill.order'
    plot_data_new <- dplyr::mutate_if(plot_data_new, is.factor, as.character)
    plot_data_new <- plot_data_new[order(plot_data_new$fill.order),]

    ## plot data
    plot <- ggplot2::ggplot(plot_data_new, ggplot2::aes(x=x, y=y, color=fill.order, fill=fill.order)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), colour = NA)

  }

  # Adjust fill colors and legend
  if(plot.type %in% c("area.edged", "area.smooth")){
    plot <- plot + ggplot2::scale_fill_manual(values = unique(plot_data[,"Col"]),labels = c(unique(plot_data[,"Fill"])))
  } else {
    if(plot.type == "line.smooth"){
      plot <- plot + ggplot2::scale_color_manual(values = unique(plot_data_new[,"Col"]),labels = c(unique(plot_data_new[,"Fill"])))
      plot <- plot + ggplot2::scale_fill_manual(values = alpha(unique(plot_data_new[,"Col"]),0.1),labels = c(unique(plot_data_new[,"Fill"])))
    } else {
      plot <- plot + ggplot2::scale_color_manual(values = unique(plot_data[,"Col"]),labels = c(unique(plot_data[,"Fill"])))
      plot <- plot + ggplot2::scale_fill_manual(values = alpha(unique(plot_data[,"Col"]),0.1),labels = c(unique(plot_data[,"Fill"])))
    }

    # legend title per color
    if(!missing(legend.names)){
      plot <-  plot + ggplot2::labs(color=legend.names)
    }
  }

  # legend title per legend.names
  if(!missing(legend.names)){
    plot <- plot + ggplot2::labs(fill=legend.names)
  }


  return(plot)
}
