style.sheet <- function(plot, xlims, ylims, colors, color.gradient, colors.opacity = 1,
                        fill.colors, fill.color.gradient, fill.colors.opacity = 1, point.shapes,
                        plot.title, plot.title.size = 16, plot.title.horizontal.position = 0.5,
                        plot.title.vertical.position = 0.5, plot.title.rotation = 0, plot.title.color = "black", plot.title.style = "plain",
                        x.axis.label, x.axis.label.size = 13, x.axis.label.horizontal.position = 0.5,
                        x.axis.label.vertical.position = 0.5, x.axis.label.rotation = 0, x.axis.label.color = "black",
                        x.axis.label.style = "plain", y.axis.label, y.axis.label.size = 13,
                        y.axis.label.horizontal.position = 0.5, y.axis.label.vertical.position = 0.5,
                        y.axis.label.rotation = 90, y.axis.label.color = "black", y.axis.label.style = "plain",
                        x.axis.text.tick.marks, x.axis.text.labels, x.axis.text.size = 12,
                        x.axis.text.horizontal.position = 0.5, x.axis.text.vertical.position = 1,
                        x.axis.text.rotation = 0, x.axis.text.color = "black", x.axis.text.style = "plain",
                        y.axis.text.tick.marks, y.axis.text.labels, y.axis.text.size = 12,
                        y.axis.text.horizontal.position = 0.5, y.axis.text.vertical.position = 1,
                        y.axis.text.rotation = 0, y.axis.text.color = "black", y.axis.text.style = "plain",
                        show.legend.by.color = TRUE, legend.title.by.color, legend.text.by.color, legend.text.by.color.order,
                        show.legend.by.fill = TRUE, legend.title.by.fill, legend.text.by.fill, legend.text.by.fill.order,
                        show.legend.by.shape = TRUE, legend.title.by.shape, legend.text.by.shape, legend.text.by.shape.order,
                        legend.position = "right", legend.columns, legend.title.size = 13, legend.title.style = "plain", legend.text.size = 12,
                        legend.text.style = "plain"){
  
  #' @title Style sheet for ggplots
  #'
  #' @description Manually adjusting graphs created with ggplot.
  #'
  #' @section Required packages:
  #' 'ggplot2'
  #'
  #' @param plot (\emph{ggplot}). A ggplot object.
  #' @param xlims (\emph{numeric/character}). An optional vector. If the x-axis-scale is \strong{numeric}, please provide two numerics specifying the range shown on the y-axis. If the data, shown on the y-axis, are \strong{characters}, please provide the names of the characters,  that shall be shown, in the desired order.
  #' @param ylims (\emph{numeric/character}). An optional vector. If the y-axis-scale is \strong{numeric}, please provide two numerics specifying the range shown on the y-axis. If the data, shown on the y-axis, are \strong{characters}, please provide the names of the characters,  that shall be shown, in the desired order.
  #' @param colors (\emph{character}). An optional vector, specifying colors for the plot.
  #' @param color.gradient (\emph{character}). An optional vector, specifying colors for a n-colour gradient.
  #' @param colors.opacity (\emph{numeric}). An optional numeric between 0-1, specifying the opacity of the colors. By default 1.
  #' @param fill.colors (\emph{character}). An optional vector, specifying colors per fill for the plot.
  #' @param fill.color.gradient (\emph{character}). An optional vector, specifying a n-colour gradient for the fill colors.
  #' @param fill.colors.opacity (\emph{numeric}). An optional numeric between 0-1, specifying the opacity of the fill colors. By default 1.
  #' @param point.shapes (\emph{numeric}). An optional vector, specifying the shapes of the points (if existing) shown in the plot.
  #' @param plot.title (\emph{character}). An optional character, defining a plot title.
  #' @param plot.title.size (\emph{numeric}). Text size of the plot title in pts. By default 16.
  #' @param plot.title.horizontal.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Horizontal position of the plot title. By default 0.5.
  #' @param plot.title.vertical.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Vertical position of the plot title. By default 1.
  #' @param plot.title.rotation (\emph{numeric}). Numeric ranging between 0-360, defining the angle of the plot title. By default 0.
  #' @param plot.title.color (\emph{character}). Color of the plot title. By default "black".
  #' @param plot.title.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the plot title. By default "plain".
  #' @param x.axis.label (\emph{character}). An optional character, adjusting x-axis label.
  #' @param x.axis.label.size (\emph{numeric}). Text size of the x-axis label in pts. By default 13.
  #' @param x.axis.label.horizontal.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Horizontal position of the x-axis label. By default 0.5.
  #' @param x.axis.label.vertical.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Vertical position of the x-axis label. By default 0.5.
  #' @param x.axis.label.rotation (\emph{numeric}). Numeric ranging between 0-360, defining the angle of the x-axis label. By default 0.
  #' @param x.axis.label.color (\emph{character}). Color of the x-axis label. By default "black".
  #' @param x.axis.label.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the x-axis label. By default "plain".
  #' @param y.axis.label (\emph{character}). An optional character, adjusting y-axis label.
  #' @param y.axis.label.size (\emph{numeric}). Text size of the y-axis label in pts. By default 13.
  #' @param y.axis.label.horizontal.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Horizontal position of the y-axis label. By default 0.5.
  #' @param y.axis.label.vertical.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Vertical position of the y-axis label. By default 0.5.
  #' @param y.axis.label.rotation (\emph{numeric}). Numeric ranging between 0-360, defining the angle of the y-axis label. By default 90.
  #' @param y.axis.label.color (\emph{character}). Color of the y-axis label. By default "black".
  #' @param y.axis.label.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the x-axis label. By default "plain".
  #' @param x.axis.text.tick.marks (\emph{numeric/character}). An optional vector, to specify which tick marks of the x-axis shall be shown.
  #' @param x.axis.text.labels (\emph{numeric/character}). An optional vector, to set the tick mark labels of the x-axis manually.
  #' @param x.axis.text.size (\emph{numeric}). Text size of the x-axis text in pts. By default 12.
  #' @param x.axis.text.horizontal.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Horizontal position of the x-axis text. By default  0.5.
  #' @param x.axis.text.vertical.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Vertical position of the x-axis text. By default 1.
  #' @param x.axis.text.rotation (\emph{numeric}). Numeric ranging between 0-360, defining the angle of the x-axis text. By default 0.
  #' @param x.axis.text.color (\emph{character}). Color of the x-axis text. By default "black".
  #' @param x.axis.text.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the x-axis text. By default "plain".
  #' @param y.axis.text.tick.marks (\emph{numeric/character}). An optional vector, to specify which tick marks of the y-axis shall be shown.
  #' @param y.axis.text.labels (\emph{numeric/character}). An optional vector, to set the tick mark labels of the y-axis manually.
  #' @param y.axis.text.size (\emph{numeric}). Text size of the y-axis text in pts. By default 12.
  #' @param y.axis.text.horizontal.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Horizontal position of the y-axis text. By default  0.5.
  #' @param y.axis.text.vertical.position (\emph{"top", "middle", "bottom", "left", "center", "right" or numeric}). Vertical position of the y-axis text. By default 1.
  #' @param y.axis.text.rotation (\emph{numeric}). Numeric ranging between 0-360, defining the angle of the y-axis text. By default 0.
  #' @param y.axis.text.color (\emph{character}). Color of the y-axis text. By default "black".
  #' @param y.axis.text.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the y-axis text. By default "plain".
  #' @param show.legend.by.color (\emph{logical}).  If \strong{TRUE}, the legend for the colors is shown, if \strong{FALSE} not. By default TRUE.
  #' @param legend.title.by.color (\emph{character}). An optional character, adjusting the title of the legend for the colors.
  #' @param legend.text.by.color (\emph{numercis/characters}). An optional vector, replacing the text of the legend for the colors.
  #' @param legend.text.by.color.order (\emph{numercis/characters}). An optional vector, changing the order of the text of the legend for the colors.
  #' @param show.legend.by.fill (\emph{logical}).  If \strong{TRUE}, the legend for the fill colors is shown, if \strong{FALSE} not. By default TRUE.
  #' @param legend.title.by.fill (\emph{character}). An optional character, adjusting the title of the legend for the fill colors.
  #' @param legend.text.by.fill (\emph{numercis/characters}). An optional vector, replacing the text of the legend for the fill colors.
  #' @param legend.text.by.fill.order (\emph{numercis/characters}). An optional vector, changing the order of the text of the legend for the fill colors.
  #' @param show.legend.by.shape (\emph{logical}).  If \strong{TRUE}, the legend for the point shapes is shown, if \strong{FALSE} not. By default TRUE.
  #' @param legend.title.by.shape (\emph{character}). An optional character, adjusting the title of the legend for the point shapes.
  #' @param legend.text.by.shape (\emph{numercis/characters}). An optional vector, replacing the text of the legend for the point shapes.
  #' @param legend.text.by.shape.order (\emph{numercis/characters}). An optional vector, changing the order of the text of the legend for the point shapes.
  #' @param legend.position (\emph{character or numeric}). The position of the legend \strong{\emph{("none", "left", "right", "bottom", "top", or two-element numeric vector)}}. By default "right".
  #' @param legend.columns (\emph{numeric}). An optional numeric defining the number of columns of the legend.
  #' @param legend.title.size (\emph{numeric}). Text size of the legend title in pts. By default 13.
  #' @param legend.title.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the legend title. By default "plain".
  #' @param legend.text.size (\emph{numeric}). Text size of the legend text in pts. By default 12.
  #' @param legend.text.style (\emph{"plain", "italic", "bold", "bold.italic"}). Text style of the legend text. By default "plain".
  #'
  #' @return
  #' @export
  #'
  #' @examples

  library(ggplot2)
  
  #############################################################################################
  ##################################### Adjusting the plot ####################################
  
  #  set plot title
  if(!missing(plot.title)){
    plot <- plot + labs(title=plot.title)
  }
  
  # x-axis label
  if(!missing(x.axis.label)){
    plot <- plot + labs(x=x.axis.label)
  }
  
  # y-axis label
  if(!missing(y.axis.label)){
    plot <- plot + labs(y=y.axis.label)
  }
  
  # legend title by color
  if(!missing(legend.title.by.color)){
    plot <- plot + labs(color=legend.title.by.color)
  }
  
  # change legend text by color, its order and colors
  if(!missing(colors)){
    if(!missing(legend.text.by.color) & !missing(legend.text.by.color.order)){
      plot <- plot + scale_color_manual(values = alpha(colors, colors.opacity), labels = legend.text.by.color, breaks = legend.text.by.color.order)
    } else if(!missing(legend.text.by.color)){
      plot <- plot + scale_color_manual(values = alpha(colors, colors.opacity), labels = legend.text.by.color)
    } else if(!missing(legend.text.by.color.order)){
      plot <- plot + scale_color_manual(values = alpha(colors, colors.opacity), breaks = legend.text.by.color.order)
    } else {
      plot <- plot + scale_color_manual(values = alpha(colors, colors.opacity))
    }
  } else if(!missing(legend.text.by.color)){
    if(!missing(legend.text.by.color.order)){
      plot <- plot + scale_color_discrete(labels = legend.text.by.color, breaks = legend.text.by.color.order)
    } else {
      plot <- plot + scale_color_discrete(labels = legend.text.by.color)
    }
  } else if(!missing(legend.text.by.color.order)){
    plot <- plot + scale_color_discrete(breaks = legend.text.by.color.order)
  }
  
  # Creates a n-colour gradients for colors and fill colors
  if(!missing(color.gradient)){
    plot <- plot + scale_color_gradientn(colors=color.gradient)
  }
  if(!missing(fill.color.gradient)){
    plot <- plot + scale_fill_gradient2(colors=fill.color.gradient)
  }
  
  # legend title per fill
  if(!missing(legend.title.by.fill)){
    plot <- plot + labs(fill=legend.title.by.fill)
  }
  
  # change legend text by fill, its order and fill colors
  if(!missing(fill.colors)){
    if(!missing(legend.text.by.fill) & !missing(legend.text.by.fill.order)){
      plot <- plot + scale_fill_manual(values = alpha(fill.colors, fill.colors.opacity), labels = legend.text.by.fill, breaks = legend.text.by.fill.order)
    } else if(!missing(legend.text.by.fill)){
      plot <- plot + scale_fill_manual(values = alpha(fill.colors, fill.colors.opacity), labels = legend.text.by.fill)
    } else if(!missing(legend.text.by.fill.order)){
      plot <- plot + scale_fill_manual(values = alpha(fill.colors, fill.colors.opacity), breaks = legend.text.by.fill.order)
    } else {
      plot <- plot + scale_fill_manual(values = alpha(fill.colors, fill.colors.opacity))
    }
  } else if(!missing(legend.text.by.fill)){
    if(!missing(legend.text.by.fill.order)){
      plot <- plot + scale_fill_discrete(values = alpha(NA, fill.colors.opacity), labels = legend.text.by.fill, breaks = legend.text.by.fill.order)
    } else {
      plot <- plot + scale_fill_discrete(labels = legend.text.by.fill)
    }
  } else if(!missing(legend.text.by.fill.order)){
    plot <- plot + scale_fill_discrete(breaks = legend.text.by.fill.order)
  }
  
  # legend title per shape
  if(!missing(legend.title.by.shape)){
    plot <- plot + labs(shape=legend.title.by.shape)
  }
  
  # change legend text by shape, its order and shapes
  if(!missing(point.shapes)){
    if(!missing(legend.text.by.shape) & !missing(legend.text.by.shape.order)){
      plot <- plot + scale_shape_manual(values = point.shapes, labels = legend.text.by.shape, breaks = legend.text.by.shape.order)
    } else if(!missing(legend.text.by.shape)){
      plot <- plot + scale_shape_manual(values = point.shapes, labels = legend.text.by.shape)
    } else if(!missing(legend.text.by.shape.order)){
      plot <- plot + scale_shape_manual(values = point.shapes, breaks = legend.text.by.shape.order)
    } else {
      plot <- plot + scale_shape_manual(values = point.shapes)
    }
  } else if(!missing(legend.text.by.shape)){
    if(!missing(legend.text.by.shape.order)){
      plot <- plot + scale_shape_discrete(labels = legend.text.by.shape, breaks = legend.text.by.shape.order)
    } else {
      plot <- plot + scale_shape_discrete(labels = legend.text.by.shape)
    }
  } else if(!missing(legend.text.by.shape.order)){
    plot <- plot + scale_shape_discrete(breaks = legend.text.by.shape.order)
  }
  
  ########### NEW: cvhange number of columns of legend
  if(!missing(legend.columns)){
    plot <- plot + guides(color = guide_legend(ncol=legend.columns),
                          fill = guide_legend(ncol=legend.columns),
                          shape = guide_legend(ncol=legend.columns))
  }
  ###########
  
  # choose if legends per colors, fill colors and point shapes are shown
  if(show.legend.by.color == FALSE){
    plot <- plot + guides(color = F)
  }
  if(show.legend.by.fill == FALSE){
    plot <- plot + guides(fill = F)
  }
  if(show.legend.by.shape == FALSE){
    plot <- plot + guides(shape = F)
  }
  
  # Set lims, tick marks and/or labels of x-axis manually
  ###### numeric
  if(is.numeric(ggplot_build(plot)$layout$panel_scales_x[[1]]$range$range)){
    if(!missing(x.axis.text.tick.marks)){
      if(!missing(x.axis.text.labels) & !missing(xlims)){
        plot <- plot + scale_x_continuous(breaks = x.axis.text.tick.marks, labels = x.axis.text.labels, limits = xlims)
      } else if(!missing(x.axis.text.labels)){
        plot <- plot + scale_x_continuous(breaks = x.axis.text.tick.marks, labels = x.axis.text.labels)
      } else if(!missing(xlims)){
        plot <- plot + scale_x_continuous(breaks = x.axis.text.tick.marks, limits = xlims)
      } else {
        plot <- plot + scale_x_continuous(breaks = x.axis.text.tick.marks)
      }
    } else if(!missing(x.axis.text.labels)){
      if(!missing(xlims)){
        plot <- plot + scale_x_continuous(labels = x.axis.text.labels, limits = xlims)
      } else {
        plot <- plot + scale_x_continuous(labels = x.axis.text.labels)
      }
    } else if(!missing(xlims)){
      plot <- plot + scale_x_continuous(limits = xlims)
    }
    ####### characters
  } else {
    if(!missing(x.axis.text.tick.marks)){
      if(!missing(x.axis.text.labels) & !missing(xlims)){
        plot <- plot + scale_x_discrete(breaks = x.axis.text.tick.marks, labels = x.axis.text.labels, limits = xlims)
      } else if(!missing(x.axis.text.labels)){
        plot <- plot + scale_x_discrete(breaks = x.axis.text.tick.marks, labels = x.axis.text.labels)
      } else if(!missing(xlims)){
        plot <- plot + scale_x_discrete(breaks = x.axis.text.tick.marks, limits = xlims)
      } else {
        plot <- plot + scale_x_discrete(breaks = x.axis.text.tick.marks)
      }
    } else if(!missing(x.axis.text.labels)){
      if(!missing(xlims)){
        plot <- plot + scale_x_discrete(labels = x.axis.text.labels, limits = xlims)
      } else {
        plot <- plot + scale_x_discrete(labels = x.axis.text.labels)
      }
    } else if(!missing(xlims)){
      plot <- plot + scale_x_discrete(limits = xlims)
    }
  }
  
  # Set lims, tick marks and/or labels of y-axis manually
  ###### numeric
  if(is.numeric(ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range)){
    if(!missing(y.axis.text.tick.marks)){
      if(!missing(y.axis.text.labels) & !missing(ylims)){
        plot <- plot + scale_y_continuous(breaks = y.axis.text.tick.marks, labels = y.axis.text.labels, limits = ylims)
      } else if(!missing(y.axis.text.labels)){
        plot <- plot + scale_y_continuous(breaks = y.axis.text.tick.marks, labels = y.axis.text.labels)
      } else if(!missing(ylims)){
        plot <- plot + scale_y_continuous(breaks = y.axis.text.tick.marks, limits = ylims)
      } else {
        plot <- plot + scale_y_continuous(breaks = y.axis.text.tick.marks)
      }
    } else if(!missing(y.axis.text.labels)){
      if(!missing(ylims)){
        plot <- plot + scale_y_continuous(labels = y.axis.text.labels, limits = ylims)
      } else {
        plot <- plot + scale_y_continuous(labels = y.axis.text.labels)
      }
    } else if(!missing(ylims)){
      plot <- plot + scale_y_continuous(limits = ylims)
    }
    ####### characters
  } else {
    if(!missing(y.axis.text.tick.marks)){
      if(!missing(y.axis.text.labels) & !missing(ylims)){
        plot <- plot + scale_y_discrete(breaks = y.axis.text.tick.marks, labels = y.axis.text.labels, limits = ylims)
      } else if(!missing(y.axis.text.labels)){
        plot <- plot + scale_y_discrete(breaks = y.axis.text.tick.marks, labels = y.axis.text.labels)
      } else if(!missing(ylims)){
        plot <- plot + scale_y_discrete(breaks = y.axis.text.tick.marks, limits = ylims)
      } else {
        plot <- plot + scale_y_discrete(breaks = y.axis.text.tick.marks)
      }
    } else if(!missing(y.axis.text.labels)){
      if(!missing(ylims)){
        plot <- plot + scale_y_discrete(labels = y.axis.text.labels, limits = ylims)
      } else {
        plot <- plot + scale_y_discrete(labels = y.axis.text.labels)
      }
    } else if(!missing(ylims)){
      plot <- plot + scale_y_discrete(limits = ylims)
    }
  }
  
  plot <- plot +
    # Title
    theme(plot.title = element_text(size = plot.title.size,
                                    hjust = plot.title.horizontal.position,
                                    vjust = plot.title.vertical.position,
                                    angle = plot.title.rotation,
                                    color = plot.title.color,
                                    face = plot.title.style),
          
          # X-axis labels
          axis.title.x = element_text(size = x.axis.label.size,
                                      hjust = x.axis.label.horizontal.position,
                                      vjust = x.axis.label.vertical.position,
                                      angle = x.axis.label.rotation,
                                      color = x.axis.label.color,
                                      face = x.axis.label.style),
          # y-axis labels
          axis.title.y = element_text(size = y.axis.label.size,
                                      hjust = y.axis.label.horizontal.position,
                                      vjust = y.axis.label.vertical.position,
                                      angle = y.axis.label.rotation,
                                      color = y.axis.label.color,
                                      face = y.axis.label.style),
          # x-axis text
          axis.text.x = element_text(size = x.axis.text.size,
                                     hjust = x.axis.text.horizontal.position,
                                     vjust = x.axis.text.vertical.position,
                                     angle = x.axis.text.rotation,
                                     color = x.axis.text.color,
                                     face = x.axis.text.style),
          # y-axis text
          axis.text.y = element_text(size = y.axis.text.size,
                                     hjust = y.axis.text.horizontal.position,
                                     vjust = y.axis.text.vertical.position,
                                     angle = y.axis.text.rotation,
                                     color = y.axis.text.color,
                                     face = y.axis.text.style),
          
          legend.position = legend.position,
          legend.title = element_text(size = legend.title.size, face = legend.title.style),
          legend.text = element_text(size = legend.text.size, face = legend.text.style))
  return(plot)
}
