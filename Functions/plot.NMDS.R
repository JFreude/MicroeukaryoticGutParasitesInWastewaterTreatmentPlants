plot.NMDS <- function(plot.data, dimensions = c(1,2), point.size = 1, group.attrib, cluster.attrib,
                             draw.ellipse.by, ellipse.size = 0.8, theme.ordination = TRUE){

  #' @title Visualization of ordinations using ggplot
  #'
  #' @description Function visualizes eigenvectors, resulting from an ordination, using ggplot. It is capable to plot the samples in different colors and shapes. Also, it can visualize groups/clusters by ellipses.
  #'
  #' @param plot.data (\emph{data.frame or matrix}). A 'data.frame' or matrix providing eigenvectors resulting from an ordination analysis with respectively one column per dimension. Dimensions must be arranged in ascending order.
  #' @param dimensions (\emph{integer}). A vector providing two integers specifying the index of two columns of plot.data that are used for the resulting plot. By default, the first and second columns/dimensions are visualized c(1,2).
  #' @param point.size (\emph{numeric}). A numeric specifying the size of all points. By default point.size = 1.
  #' @param group.attrib (\emph{characters/numercis}). An optional vector grouping the sample. The groups are visualized by colors. The length of the vector must match the number of rows of \strong{\emph{plot.data}}.
  #' @param cluster.attrib (\emph{characters/numercis}).An optional vector clustering the sample. The clusters are visualized by shapes. The length of the vector must match the number of rows of \strong{\emph{plot.data}}.
  #' @param draw.ellipse.by (\emph{"Group", "Cluster"}). Optional visualization of polygons encapsulating the groups or the clusters using ggplot's stat_ellipse function. Consequentially, the argument can only be set to \strong{"Group"} if the argument \strong{\emph{group.attrib}} is defined. It can only be set to \strong{"Cluster"} if the argument \strong{\emph{cluster.attrib}} is defined.
  #' @param ellipse.size (\emph{numeric}). A numeric ranging between 0-1 (0 < value < 1) specifying the size of the ellipses. By default ellipse.size = 0.8.
  #' @param theme.ordination (\emph{logical}). If \strong{TRUE}, the theme is adjusted, only the outer lines and a horizontal and vertical zero line are shown, if \strong{FALSE} the theme is not adjusted. By default theme.ordination = TRUE.
  #'
  #' @section Packages the function bases on
  #' 'ggplot2'

  #############################################################################################
  #################################### Check prerequisites ####################################

  #----------------- Check if all needed functions and packages are available ----------------#

  if(!"ggplot2" %in% rownames(installed.packages())){
    stop("Please install the package 'ggplot2'.")
  }

  #--------------------------- Check if needed arguments are valid ---------------------------#

  # Check if 'plot.data' are provided
  if(missing(plot.data) | !(class(plot.data) == "data.frame" | class(plot.data) == "matrix")){
    stop("Please provide a 'data.frame' or 'matrix' with eigenvectors/dimensions resulting from an ordination analysis. Dimensions must be arranged in ascending order, one per column.")
  } else {
    plot.data <- data.frame(plot.data)
  }

  # Check if the values of the argument 'dimensions' are valid
  if(length(dimensions) != 2 | class(dimensions) != "numeric"){
    stop("Please provide two integers, specifying the columns/dimensions that shall be plotted.")
  } else if(any(dimensions%%1 != 0)){
    stop("Please provide two integers, specifying the columns/dimensions that shall be plotted.")
  } else if(any(dimensions > ncol(plot.data))){
    stop("One or both of the provided column/dimension indices are higher than the number of columns of the plot.data.")
  }

  # Check if the value of the argument 'point.size' is valid
  if(length(point.size) != 1 | class(point.size) != "numeric"){
    stop("Please provide one numeric, specifying the point size.")
  }

  #-------------------------- Check if optional arguments are valid --------------------------#
  # Check if the value of the argument 'group.attrib' is valid
  if(missing(group.attrib)){
    group.attrib <- NULL
  } else if(length(group.attrib) != nrow(plot.data)){
    stop("The length of the vector 'group.attrib' must match the number of rows of 'plot.data'")
  }

  # Check if the value of the argument 'cluster.attrib' is valid
  if(missing(cluster.attrib)){
    cluster.attrib <- NULL
  } else if(length(cluster.attrib) != nrow(plot.data)){
    stop("The length of the vector 'cluster.attrib' must match the number of rows of 'plot.data'")
  }

  # Check if the value of the argument 'draw.ellipse.by' is valid
  if(missing(draw.ellipse.by)){
    draw.ellipse.by <- NULL
  } else if(!draw.ellipse.by %in% c("Group", "Cluster")){
    stop("Wrong value for the argument 'draw.ellipse.by' is provided. Please choose between 'Group' and 'Cluster'.")
  } else if(draw.ellipse.by == "Group" & is.null(group.attrib)){
    warning("Ellipses by 'Group' can not be drawn, since the argument 'group.attrib' is not defined. No ellipses are shown in the plot.")
    draw.ellipse.by <- NULL
  } else if(draw.ellipse.by == "Cluster" & is.null(cluster.attrib)){
    warning("Ellipses by 'Cluster' can not be drawn, since the argument 'cluster.attrib' is not defined. No ellipses are shown in the plot.")
    draw.ellipse.by <- NULL
  }

  # Check if the value of the argument 'ellipse.size' is valid
  if(!is.numeric(ellipse.size) | length(ellipse.size) != 1){
    stop("Please provide one numeric ranging between 0-1 to specify the size of the shown ellipses.")
  } else if(!(ellipse.size >= 0 & ellipse.size <= 1)){
    stop("The length of the vector 'cluster.attrib' must match the number of rows of 'plot.data'")
  }

  # Check if the value of the argument 'theme.ordination' is valid
  if(!theme.ordination %in% c(TRUE, FALSE)){
    stop("Wrong value for the argument 'theme.ordination' is provided. Please choose between TRUE and FALSE.")
  }

  #############################################################################################
  ############################### Create dataframe for plotting ###############################

  # Create dataframe with group and/or cluster (if given)
  if(is.null(group.attrib) & is.null(cluster.attrib)){ # No Group or Cluster
    plot.coords <- data.frame(plot.data[,dimensions], Group="", Cluster="")
    plot.version <- "basic"
  } else if(!is.null(group.attrib) & !is.null(cluster.attrib)){ # Group and Cluster
    plot.coords <- data.frame(plot.data[,dimensions], Group=group.attrib, Cluster=cluster.attrib)
    plot.version <- "groups_clusters"
  } else if(!is.null(group.attrib)){ # Group
    plot.coords <- data.frame(plot.data[,dimensions], Group=group.attrib, Cluster="")
    plot.version <- "groups"
  } else if(!is.null(cluster.attrib)){ # Cluster
    plot.coords <- data.frame(plot.data[,dimensions], Group="", Cluster=cluster.attrib)
    plot.version <- "clusters"
  }

  #############################################################################################
  ########################################### Plots ###########################################

  #---------------------------------------- Basic plot ---------------------------------------#
  plot <- ggplot2::ggplot(data=plot.coords, ggplot2::aes(x=plot.coords[,1], y=plot.coords[,2], color=Group))

  #-------------------------------------- Adding ellipse -------------------------------------#
  if(!is.null(draw.ellipse.by)){
    if(draw.ellipse.by == "Group"){
      plot <- plot + ggplot2::stat_ellipse(geom = "polygon", linetype = "blank", ggplot2::aes(fill = Group),
                                  level = ellipse.size, show.legend = FALSE)
    }else if(draw.ellipse.by == "Cluster"){
      plot <- plot + ggplot2::stat_ellipse(geom = "polygon", linetype = "blank", ggplot2::aes(fill = Cluster),
                                  level = ellipse.size, show.legend = FALSE)
    }
  }

  #--------------------------------------- Adding point --------------------------------------#

  plot <- plot + ggplot2::geom_point(ggplot2::aes(shape=Cluster), size = point.size)

  # Shape points manually if there are more than 6 shapes
  if(length(unique(plot.coords$Cluster)) > 6){
    plot <- plot + ggplot2::scale_shape_manual(values=seq(unique(plot.coords$Cluster)))
  }

  #-------------------------------------- Adjusting axis -------------------------------------#

  # Rename axis
  plot <- plot + ggplot2::labs(x=names(plot.coords)[1], y=names(plot.coords)[2])

  # Adjust legends, Group and Clusters are always given to adjust colors easy but if the columns are empty, also no legend should be shown
  if(plot.version == "basic"){ # No group or cluster
    plot <- plot + ggplot2::guides(color=F, shape=F)
  } else if (plot.version == "groups"){
    plot <- plot + ggplot2::guides(shape=F)
  } else if (plot.version == "clusters"){
    plot <- plot + ggplot2::guides(color=F)
  }

  #--------------------------------------- Adjust theme --------------------------------------#
  if(theme.ordination){
    plot <- plot + ggplot2::geom_hline(yintercept = 0, colour="black", linetype="dotted") +
      ggplot2::geom_vline(xintercept = 0, colour="black", linetype="dotted") +
      ggplot2::theme_linedraw()+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank())
  }

  #-------------------------------------- Set axis range -------------------------------------#

  if(is.null(draw.ellipse.by)){
    xmin <- min(ggplot2::ggplot_build(plot)$data[[1]]$x)
    xmax <- max(ggplot2::ggplot_build(plot)$data[[1]]$x)
    ymin <- min(ggplot2::ggplot_build(plot)$data[[1]]$y)
    ymax <- max(ggplot2::ggplot_build(plot)$data[[1]]$y)
  } else {
    xmin <- min(c(ggplot2::ggplot_build(plot)$data[[1]]$x, ggplot2::ggplot_build(plot)$data[[2]]$x))
    xmax <- max(c(ggplot2::ggplot_build(plot)$data[[1]]$x, ggplot2::ggplot_build(plot)$data[[2]]$x))
    ymin <- min(c(ggplot2::ggplot_build(plot)$data[[1]]$y, ggplot2::ggplot_build(plot)$data[[2]]$y))
    ymax <- max(c(ggplot2::ggplot_build(plot)$data[[1]]$y, ggplot2::ggplot_build(plot)$data[[2]]$y))
  }

  plot <- plot + ggplot2::xlim(xmin - (xmax - xmin)*0.02, xmax + (xmax - xmin)*0.02) +
    ggplot2::ylim(ymin - (ymax - ymin)*0.02, ymax + (ymax - ymin)*0.02)

  return(plot)
}
