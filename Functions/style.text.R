style.text <- function(names, italic.text, p.values){
  
  #' @title Italic vs. plain text style for taxa/metadata names plotted with 'geom_text_repel'
  #'
  #' @description The function is able to format taxa/metadata names in a way, that they are plotted in italic or plain while visualization with ggplot's function 'geom_text_repel'. Further, significant levels can be visualized by stars.
  #'
  #' @param names (\emph{character}). A vector providing taxa/metadata names as characters.
  #' @param italic.text (\emph{logical}). A vector providing logical values (TRUE/FALSE) for each taxon/metadata factor of \strong{\emph{'names'}}. If \strong{TRUE}, taxon/factor is plotted in italic text style if \strong{FALSE}, taxon/factor is plotted in plain text style.
  #' @param p.values (\emph{numeric}). An optional vector providing p-values for each taxon/metadata factor of \strong{\emph{'names'}}. Significant levels are visualized by stars. If p-value <= 0.001, 3 stars, if p-value <= 0.01 2 stars, if p-value < 0.05  1 star and if p-value >= 0.05 zero stars are plotted.
  

  #############################################################################################
  #################################### Check prerequisites ####################################

  #--------------------------- Check if needed arguments are valid ---------------------------#

  # Check if 'names' are provided
  if(missing(names) | class(names) != "character"){
    stop("Please provide a vector with names as characters.")
  }

  # Check if the argument 'italic.text' is valid
  if(missing(italic.text) | class(italic.text) != "logical" | length(italic.text) != length(names)){
    stop("Please provide a vector with logicals (italic.text = TRUE/FALSE) specifying the text style for each taxon/metadata factor.")
  }

  # Check if the argument 'p.values' is valid
  if(missing(p.values)){
    p.values <- NULL
  } else if(class(p.values) != "numeric" | length(p.values) != length(names)){
    stop("If desired, please provide a vector with numerics specifying the p-values for each taxon/metadata factor.")
  }

  #############################################################################################
  ######################################## Format text ########################################

  # Open list for results
  formatted.names <- list()

  # Iterates over each name
  for(name.id in seq(names)){
    name <- names[name.id]

    #--------------------------------------- No p-values ---------------------------------------#
    if(is.null(p.values)){
      if(italic.text[name.id]){# italic
        formatted.names[[name.id]] <- bquote(italic(.(name)))
      } else {# plain
        formatted.names[[name.id]] <- bquote(plain(.(name)))
      }
    #----------------------------------------- p-values ----------------------------------------#
    } else {
      pval <- p.values[name.id]
      if(italic.text[name.id]){# italic
        if (pval <= 0.001 ){
          formatted.names[[name.id]] <- bquote(italic(.(name)) ~ "***")
        } else if (pval <= 0.01){
          formatted.names[[name.id]] <- bquote(italic(.(name)) ~ "**")
        } else if (pval < 0.05){
          formatted.names[[name.id]] <- bquote(italic(.(name)) ~ "*")
        } else {
          formatted.names[[name.id]] <- bquote(italic(.(name)))
        }
      } else {
        if (pval <= 0.001 ){# plain
          formatted.names[[name.id]] <- bquote(plain(.(name)) ~ "***")
        } else if (pval <= 0.01){
          formatted.names[[name.id]] <- bquote(plain(.(name)) ~ "**")
        } else if ((pval < 0.05)&(pval > 0.01)){
          formatted.names[[name.id]] <- bquote(plain(.(name)) ~ "*")
        } else {
          formatted.names[[name.id]] <- bquote(plain(.(name)))
        }
      }
    }
  }

  return(unlist(formatted.names))
}
