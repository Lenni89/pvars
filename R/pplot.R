

#### Generic functions and their S3 methods ####
# CRAN manual: https://cran.r-project.org/doc/manuals/R-exts.html#Generic-functions-and-methods
# Roxygen: https://r-pkgs.org/man.html#man-s3

#' @title Coerce into a '\code{pplot}' object 
#' @description Coerce into a '\code{pplot}' object. This function allows (1) to
#'   overlay and colorize multiple plots of IRF and confidence bands in a single 
#'   '\code{ggplot}' object and (2) to overwrite shock and variable names in
#'   verbatim LaTeX math mode suitable for the export via \strong{tikzDevice}. 
#' @details \code{\link{as.pplot}} is used as an intermediary in the '\code{pplot}' 
#'   functions to achieve compatibility with different classes. Equivalently to 
#'   \code{\link{as.pvarx}} for panels of \eqn{N} VAR objects, \code{\link{as.pplot}} 
#'   summarizes panels of \eqn{G} plot objects that have been returned 
#'   from the '\code{plot}' method for class '\code{svarirf}' or '\code{sboot}'. 
#'   If the user wishes to extend the compatibility of the '\code{pplot}' 
#'   functions with further classes, she may simply specify accordant 
#'   \code{\link{as.pplot}}-methods instead of altering the original 
#'   '\code{pplot}' functions.
#'   
#' @param ... A single ggplot or list(s) of ggplots to be transformed.  
#' @param names_k Vector. Names of the variables \eqn{k=1,\ldots,K}. 
#'   If \code{NULL} (the default), the names of the first plot are reused. 
#'   For LaTeX exports, use e.g. \code{paste0("y_{ ", 1:dim_K, " }")}.
#' @param names_s Vector. Names of the shocks \eqn{s=1,\ldots,S}. If \code{NULL} 
#'   (the default), the names of the first plot are reused. 
#'   For LaTeX exports, use e.g. \code{paste0("\\\epsilon_{ ", 1:dim_S, " }")}.
#' @param names_g Vector. Names of the layered plots \eqn{g=1,\ldots,G}. If \code{NULL} 
#'   (the default), the names of the plots given in \code{...} are reused. 
#' @param color_g Vector. Colors of the layered plots \eqn{g=1,\ldots,G}. 
#' @param shape_g Vector. Shapes of the layered plots \eqn{g=1,\ldots,G}, 
#'   see \strong{ggplot2}'s \code{\link[ggplot2]{geom_point}} for shapes. 
#'   If \code{NULL} (the default), no points will be set on the IRF-lines.
#' @param width_L Integer. The \code{\link[ggplot2:aes_linetype_size_shape]{linewidth}} 
#'   of the lines tracing the ribbon edges. Note that each line is centered around 
#'   the ribbon edge, thus exceeding the interval endpoints by \code{width_L}\eqn{/2}.
#'   If \code{NULL} (the default), no lines will be set on the edges.
#' @param color_L Character. The color of the lines tracing the edges. If \code{NULL} 
#'   (the default), the plot-specific \code{color_g} will be used instead.
#' @param n.rows Integer. Number of rows in \code{\link[ggplot2]{facet_wrap}}. 
#'   If \code{NULL} (the default), the dimensions of the facet plots given 
#'   in \code{...} are reused.
#' @param scales Character. Should scales be fixed (\code{"fixed"}), 
#'  free (\code{"free"}), or free in one dimension (\code{"free_x"}, 
#'  \code{"free_y"}, the default)? See \code{\link[ggplot2]{facet_wrap}}.
#' @param Latex Logical. If \code{TRUE}, the arrows of the facet labels are 
#'   written in LaTeX math mode.
#' 
#' @return A list of class '\code{pplot}'. Objects of this class contain the elements:
#'   \item{F.plot}{'\code{ggplot}' object for the merged plot.}
#'   \item{L.plot}{List of '\code{ggplot}' objects containing all \eqn{G} plots.}
#'   \item{args_pplot}{List of characters and integers indicating the
#'                     specifications used for creating \code{F.plot}.}
#' @seealso \code{\link{PP}}, \code{\link{irf.pvarx}}, \code{\link{pid.dc}}, and 
#'   \code{\link{id.iv}} for further examples of edited plots, in particular for 
#'   subset and reordered facet plots with reshaped facet dimensions.
#' @example inst/examples/pplot.R
#' @export
#' 
as.pplot <- function(..., names_k=NULL, names_s=NULL, names_g=NULL, 
  color_g=NULL, shape_g=NULL, width_L=NULL, color_L=NULL, 
  n.rows=NULL, scales="free_y", Latex=FALSE) UseMethod("as.pplot")


#' @method as.pplot default
#' @importFrom scales hue_pal
#' @export
as.pplot.default <- function(..., names_k=NULL, names_s=NULL, names_g=NULL, 
  color_g=NULL, shape_g=NULL, width_L=NULL, color_L=NULL, 
  n.rows=NULL, scales="free_y", Latex=FALSE){
  # defaults if overwritten as NULL
  if( is.null(scales) ){ scales = "free_y" }
  if( is.null(Latex)  ){ Latex  = FALSE }
  
  # try to homogenize all listed elements to ggplot and ggplot_built objects
  aux_ggunlist <- function(x){
    L.out = list()
    for(j in 1:length(x)){
      if(inherits(x[[j]], "ggplot")){
        L.out = c(L.out, x[j])
      }else if(inherits(x[[j]], "pplot")){
        L.out = c(L.out, x[[j]]$L.plot)
      }else if(is.list(x[[j]])){
        L.out = c(L.out, aux_ggunlist(x[[j]]))
      }else{
        stop("Arguments are not objects of suitable class!")
      }
    }
    return(L.out)
  }
  
  x      = list(...)
  L.plot = aux_ggunlist(x)
  dim_G  = length(L.plot)  # number of plotted groups
  if( !all(sapply(L.plot, FUN=function(x) inherits(x, "ggplot"))) ){
     stop("Arguments are not objects of suitable class!") }
  L.built = lapply(L.plot, FUN=function(x) ggplot2::ggplot_build(x))
  
  # names for variables, for shocks, and for the header of each IRF
  df_lay = L.built[[1]]$layout$layout  # get names from first ggplot object
  R.labs = do.call("rbind", sapply(df_lay$variable, FUN=function(x)
    strsplit(as.character(x), split=c(" %->% ", " \\rightarrow "))))
  
  if( is.null(names_k) ){ 
    names_k = unique(R.labs[, 2])
    df_lay$names_k = R.labs[, 2]
  }else{
    df_lay$names_k = names_k[df_lay$ROW]
  }
  
  if( is.null(names_s) ){ 
    names_s = unique(R.labs[, 1])
    df_lay$names_s = R.labs[, 1]
  }else{
    df_lay$names_s = names_s[df_lay$COL]
  }
  
  if(Latex){ 
    names_IRF = paste0("$ ", df_lay$names_s, " \\rightarrow ", df_lay$names_k, " $")
    names(names_IRF) = names_IRF
    label_IRF = as_labeller(names_IRF, default=label_value)
  }else{ 
    names_IRF = paste0(df_lay$names_s, " %->% ", df_lay$names_k)
    names(names_IRF) = names_IRF
    label_IRF = as_labeller(names_IRF, default=label_parsed)
  }
  
  # names for groups, their colors, and shapes
  if( is.null(names_g) ){ names_g = names(L.plot) }else{ names(L.plot) = names_g }
  if( dim_G != length(unique(names_g)) ){
    stop("The layered plots do not have unique names. Please provide ", dim_G, 
          " unique names via the attribute 'names_g' or label all input plots.") 
  }
  n.shps = dim_G-length(shape_g)
  n.clrs = dim_G-length(color_g)
  R.shps = c(shape_g, if(n.shps > 0) rep(0:6, length.out=n.shps))
  R.clrs = c(color_g, if(n.clrs > 0) scales::hue_pal()(n=n.clrs))  # ggplot default colors
  R.rgb  = grDevices::col2rgb(R.clrs)/3 + 170  # 'darkgray' = 'black'/3 + 169 
  R.fill = grDevices::rgb(R.rgb[1, ], R.rgb[2, ], R.rgb[3, ], maxColorValue=255)
  names(names_g) = names_g
  names(R.shps)  = names_g
  names(R.clrs)  = names_g
  names(R.fill)  = names_g
  names_g = factor(names_g, levels=names_g)  # ordering of names_g and of layers must be identical 
  
  # gather data from plots in data.frame
  L.irf = L.cbs = list()
  for(g in 1:dim_G){
    n.layers = length(L.built[[g]]$data)
    if(n.layers==2){
      # from plot.svarirf() for IRF or PP
      L.irf[[g]] = L.built[[g]]$data[[1]]
      L.cbs[[g]] = NULL
      L.irf[[g]]$variable = names_IRF[match(L.irf[[g]]$PANEL, df_lay$PANEL)]
      L.irf[[g]]$plot_g = names_g[g]
    }else{
      # from plot.sboot() for IRF with confidence bands
      L.irf[[g]] = L.built[[g]]$data[[2]]
      L.cbs[[g]] = L.built[[g]]$data[[1]]
      L.irf[[g]]$variable = names_IRF[match(L.irf[[g]]$PANEL, df_lay$PANEL)]
      L.cbs[[g]]$variable = names_IRF[match(L.cbs[[g]]$PANEL, df_lay$PANEL)]
      L.irf[[g]]$plot_g = names_g[g]
      L.cbs[[g]]$plot_g = names_g[g]
    }
  }
  df_irf = do.call("rbind", L.irf)
  df_cbs = do.call("rbind", L.cbs)
  
  # stfu R CMD check vs. ggplot2 (common practice, aes_ is deprecated)
  y = ymin = ymax = alpha = group = plot_g = NULL  # 'group' identifies different confidence bands
  
  # plot IRF with optional confidence bands
  R.guide = guide_legend(override.aes=list(shape=R.shps, fill=R.fill))  # aesthetics to merge in color legend
  F.plot  = ggplot() +
    ## geoms for df_cbs ##
    {if(!is.null(df_cbs)) 
      geom_ribbon(data = df_cbs, 
                  aes(x=x, ymin=ymin, ymax=ymax, alpha=alpha, 
                  fill=plot_g, group=interaction(group, plot_g)))} +
    {if(!is.null(df_cbs) & !is.null(width_L)) if(!is.null(color_L))
      geom_line(data = df_cbs, 
                aes(x=x, y=ymin, group=interaction(group, plot_g)), 
                show.legend=FALSE, linewidth=width_L, color=color_L) else 
      geom_line(data = df_cbs, 
                aes(x=x, y=ymin, color=plot_g, group=interaction(group, plot_g)), 
                show.legend=FALSE, linewidth=width_L)} +
    {if(!is.null(df_cbs) & !is.null(width_L)) if(!is.null(color_L))
      geom_line(data = df_cbs, 
                aes(x=x, y=ymax, group=interaction(group, plot_g)), 
                show.legend=FALSE, linewidth=width_L, color=color_L) else 
      geom_line(data = df_cbs, 
                aes(x=x, y=ymax, color=plot_g, group=interaction(group, plot_g)), 
                show.legend=FALSE, linewidth=width_L)} +
    ## geoms for df_irf ###
    geom_line(data = df_irf, aes(x=x, y=y, color=plot_g)) +
    {if(!is.null(shape_g)) 
      geom_point(data = df_irf, aes(x=x, y=y, color=plot_g, shape=plot_g))} +
    geom_hline(yintercept=0, color="red") + 
    facet_wrap(~factor(variable, levels=names_IRF), nrow=n.rows, scales=scales, labeller=label_IRF) + 
    ## scales ##
    {if(!is.null(df_cbs)) 
      scale_alpha_identity()} + 
    {if(!is.null(df_cbs)) 
      scale_fill_manual( labels=names(names_g), values=R.fill)} +
    scale_color_manual(  labels=names(names_g), values=R.clrs)  +
    {if(!is.null(shape_g)) 
      scale_shape_manual(labels=names(names_g), values=R.shps)} +
    ## layout ##
    guides(color=R.guide, shape="none", fill="none", alpha="none") +
    labs(x="Horizon", y="Response", color="Group") +
    theme_bw()
  
  # return result
  if( is.null(shape_g) ){ R.shps = NULL }
  args_pplot = list(names_k=names_k, names_s=names_s, names_g=names_g, 
    color_g=R.clrs, shape_g=R.shps, width_L=width_L, color_L=color_L,
    n.rows=n.rows, scales=scales, Latex=Latex)
  result = list(F.plot=F.plot, L.plot=L.plot, args_pplot=args_pplot)
  class(result) = "pplot"
  return(result)
}


#' @method as.pplot pplot
#' @export
as.pplot.pplot <- function(..., names_k=NULL, names_s=NULL, names_g=NULL, 
  color_g=NULL, shape_g=NULL, width_L=NULL, color_L=NULL,
  n.rows=NULL, scales=NULL, Latex=NULL){
  # gather arguments
  x = list(...)[[1]]
  args_pplot = list(names_k=names_k, names_s=names_s, names_g=names_g, 
       color_g=color_g, shape_g=shape_g, width_L=width_L, color_L=color_L,
       n.rows=n.rows, scales=scales, Latex=Latex)
  is_NULL = c(is.null(names_k), is.null(names_s), is.null(names_g), 
       is.null(color_g), is.null(shape_g), is.null(width_L), is.null(color_L),
       is.null(n.rows), is.null(scales), is.null(Latex))
  is_same = sapply(1:length(args_pplot), FUN=function(j) identical(args_pplot[j], x$args_pplot[j]))
  
  # compare arguments and return
  if( all(is_NULL | is_same) ){
    return(x)
    
  }else{
    # overwrite with given arguments
    if( is.null(names_k) ){ names_k = x$args_pplot$names_k }
    if( is.null(names_s) ){ names_s = x$args_pplot$names_s }
    if( is.null(names_g) ){ names_g = x$args_pplot$names_g }
    if( is.null(color_g) ){ color_g = x$args_pplot$color_g }
    if( is.null(shape_g) ){ shape_g = x$args_pplot$shape_g }
    if( is.null(width_L) ){ width_L = x$args_pplot$width_L }
    if( is.null(color_L) ){ color_L = x$args_pplot$color_L }
    if( is.null(n.rows)  ){ n.rows  = x$args_pplot$n.rows  }
    if( is.null(scales)  ){ scales  = x$args_pplot$scales  }
    if( is.null(Latex)   ){ Latex   = x$args_pplot$Latex   }
    
    return(as.pplot.default(..., names_k=names_k, names_s=names_s, names_g=names_g, 
      color_g=color_g, shape_g=shape_g, width_L=width_L, color_L=color_L,
      n.rows=n.rows, scales=scales, Latex=Latex))
  }
}



#### S3 methods for objects of class 'pplot' ####

#' @method plot pplot
#' @export
#' 
plot.pplot <- function(x, ...){ plot(x$F.plot) }


#' @method print pplot
#' @export
#' 
print.pplot <- function(x, ...){ print(x$F.plot) }


