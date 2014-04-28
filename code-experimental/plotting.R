#' @rdname assorted-methods
#' @aliases geta4aLatticeOptions
#' @section geta4aLatticeOptions:
#' Set lattice options for a4a plots.
#' @examples
#' # set up grapical options
#'  lattice.options(default.theme = geta4aLatticeOptions())
geta4aLatticeOptions <- function(...) 
{
  opts <- standard.theme(color = FALSE)
  opts <- modifyList(opts, list(strip.background = list(col = grey(0.95))))
  opts <- modifyList(opts, list(plot.symbol = list(pch = 16, col = "black", cex = 0.7)))
  opts <- modifyList(opts, list(axis.line = list(col = grey(0.8))))
  opts <- modifyList(opts, list(strip.border = list(col = grey(0.8))))

  opts
}

#' adds bubbles to a plot - haha!
#'
#'
#' @param x the x location of the bubbles
#' @param y the y location of the bubbles
#' @param v a vector of bubble sizes
#' @param scale a scaling factor 
#' @param ...  extra arguments to plot and points
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{extractData} is intended to be used internally
#' @author anders nielson \email{an@@dtu.aqua.dk}
bp <- function(x, y, v, scale = 3, ...)
{
  ref <- qnorm(c(.75,.975))
  col <- c(colorRampPalette(c("blue","lightblue"))(3), colorRampPalette(c("tomato2","red"))(3))
  vref <- as.numeric(factor(sign(v) * ifelse(abs(v) < ref[2], ifelse(abs(v) < ref[1], 1, 2), 3)))
  plot(x, y, cex = sqrt(abs(v)) * scale, 
       #bg = ifelse(v>0,'tomato2','lightblue'), 
       col = col[vref], 
       pch = 16, ...)
  #vmin <- order(v)[1:2]
  #vmax <- rev(order(v))[1:2]  
  #points(x[vmin], y[vmin], col = 'white', pch="-")
  #points(x[vmax], y[vmax], col='white', pch="+")

  if (any(v < -2)) text(x[v < -2], y[v < -2], "< -2",  col = 'white', cex = 0.85)
  if (any(v > 2)) text(x[v > 2], y[v > 2], "> +2", col='white', cex = 0.85)

  if (any(v < -.68 & v > -2)) points(x[v < -.68 & v > -2], y[v < -.68 & v > -2], col = 'blue', pch="-")
  if (any(v > .68 & v < 2)) points(x[v > .68 & v < 2], y[v > .68 & v < 2], col='red', pch="+")

  #if (!missing(ref)) sapply(ref, function(z) points(x, y, cex = sqrt(abs(z)) * scale, col = "white"))
  invisible(NULL)
}


#' Extracts the catch and survey data, makes useful summaries
#' and places them in an environment
#'
#'
#' @param stock an FLStock object containing catch and stock information
#' @param index an FLIndex object containing survey indices 
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{extractData} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
matrix.plot <- 
function (x, 
          xlim = 0.5 + c(0, di[2]), 
          ylim = 0.5 + c(di[1], 0),
          xlab = "", ylab = "",
          cols = c("red","blue"), text = FALSE, xmin = 1, ymin = 1, ...) 
{

  di <- dim(x)

  xlim <- xlim + xmin - 1
  ylim <- ylim + ymin - 1

  df <- data.frame(i = c(row(x)), j = c(col(x)), x = c(x))
  df $ i <- df $ i + ymin - 1
  df $ j <- df $ j + xmin - 1  

  
  xx <- df $ x
  rx <- range(xx, finite = TRUE)
  nn <- 100
  n0 <- min(nn, max(0, round((0 - rx[1])/(rx[2] - rx[1]) * nn)))
  levelplot(
    df$x ~ df$j + df$i,
    sub = "", xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
    aspect = "iso", colorkey = FALSE, 
    col.regions = colorRampPalette(rev(cols))(16), cuts = 15, 
    par.settings = list(background = list(col = "transparent")), 
    panel = function(x, y, z, subscripts, at, ..., col.regions) 
      {
        x <- as.numeric(x)
        y <- as.numeric(y)
        numcol <- length(at) - 1
        num.r <- length(col.regions)
        col.regions <- if (num.r <= numcol) {
          rep(col.regions, length = numcol)
        } else {
          col.regions[1 + ((1:numcol - 1) * (num.r - 1))%/%(numcol - 1)]
        }
        zcol <- rep.int(NA_integer_, length(z))
        for (i in seq_along(col.regions)) {
          zcol[!is.na(x) &  !is.na(y) & !is.na(z) & at[i] <= z & z < at[i + 1]] <- i
        }
        wh <- grid::current.viewport()[c("width", "height")]
        wh <- c(grid::convertWidth(wh$width, "inches", valueOnly = TRUE), grid::convertHeight(wh$height, "inches", valueOnly = TRUE)) * par("cra")/par("cin")
        pSize <- wh/di
        pA <- prod(pSize)
        p1 <- min(pSize)
        lwd <- if (p1 < 2 || pA < 6)  0.01 else if (p1 >= 4) 1 else if (p1 > 3)  0.5 else 0.2
        grid.rect(x = x, y = y, width = 1, height = 1, default.units = "native", gp = gpar(fill = col.regions[zcol],  lwd = lwd, col = if (lwd < 0.01) NA))
        # plot value inside square, scale size depending on maximum dimension
        if (text) grid.text(sprintf("%.2f", z), x = x, y = y, default.units = "native", gp = gpar(col = "darkblue", cex = 0.8 * (10 / max(di))^.4 ))
     }, ...
  )
}

#' Extracts the catch and survey data, makes useful summaries
#' and places them in an environment
#'
#'
#' @param stock an FLStock object containing catch and stock information
#' @param index an FLIndex object containing survey indices 
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{extractData} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
  panel.3d.levelplot <-
  function(x, y, z, rot.mat, distance, zlim.scaled, at, drape = TRUE, shade = FALSE, cols, ...) 
  {
    zrng <- 0.001 * diff(zlim.scaled) 
    ## vertical range of 2nd surface (>0)     
    z.scaled <- (z - min(z))/diff(range(z))     
    at.scaled <- (at - min(z))/diff(range(z)) 
    new.level <- zlim.scaled[2]    
    new.z <- new.level + zrng * (z.scaled)     
    new.at <- new.level + zrng * (at.scaled)     
    panel.3dwire(x, y, new.z, at = new.at,
                 col = grey(.9),
                 rot.mat = rot.mat, distance = distance,
                 shade = FALSE,
                 drape = TRUE,
                 zlim.scaled = zlim.scaled,
                 alpha = 1, ...)
    dots <- list(...)
    dots $ col.regions <- paste0(colorRampPalette(rev(cols))(100), "AA")
    dots $ background <- "transparent"
    args <- c(dots, list(x=x, y=y, z=z, rot.mat = rot.mat, distance = distance,
                 zlim.scaled = zlim.scaled, at = at, drape = FALSE, shade = FALSE))
    do.call(panel.3dwire, args)
  }

#' Extracts the catch and survey data, makes useful summaries
#' and places them in an environment
#'
#'
#' @param stock an FLStock object containing catch and stock information
#' @param index an FLIndex object containing survey indices 
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{extractData} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
  plotError <- function(x, y, y.se, FUN = function(x) x, ylab="", xlab="", add = FALSE, cols) {
    ylim <- c(0, FUN(max(y + 2*y.se)))
    if(is.infinite(ylim[2])) ylim <- c(0, max(y + 2 * y.se))
    if (!add) plot(x, FUN(y), ylab=ylab, xlab=xlab, type='n', ylim = ylim, las=1) 
    polygon(c(x, rev(x)), FUN(c(y - 2.58*y.se, rev(y + 2.58*y.se))), border=FALSE, col=cols[1], density=NA) 
    polygon(c(x, rev(x)), FUN(c(y - 1.96*y.se, rev(y + 1.96*y.se))), border=FALSE, col=cols[2], density=NA)
    polygon(c(x, rev(x)), FUN(c(y - .675*y.se, rev(y + .675*y.se))), border=FALSE, col=cols[3], density=NA)
    lines(x, FUN(y), lwd=2)
  }

