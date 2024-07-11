# auxiliary functions for drawing ellipses in igraph as node shapes

xyAngle <- function(x,
                    y=NULL,
                    directed=FALSE,
                    deg=TRUE,
                    origin.x=0,
                    origin.y=0,
                    ...){
  ## Get angle from zero to given x,y coordinates
  if (length(y) == 0) {
    y <- x[,2];
    x <- x[,1];
  }
  if (length(ncol(origin.x)) > 0) {
    origin.y <- origin.x[,2];
    origin.x <- origin.x[,1];
  }
  out <- base::atan2(y - origin.y,
                     x - origin.x);
  if (!directed) {
    out <- out %% pi;
  }
  if (deg) {
    out <- out * 180 / pi;
  }
  out;
}


#' Draw ellipse
#'
#' Draw ellipse
#'
#' This function draws an ellipse centered on the given coordinates,
#' rotated the given degrees relative to the center point, with give
#' x- and y-axis radius values.
#'
#' @return invisible list of x,y coordinates
#'
#' @param x,y `numeric` coordinates, where x can be a two-column numeric
#'    matrix of x,y coordinates.
#' @param a,b `numeric` values indicating x- and y-axis radius, before
#'    rotation if `angle` is non-zero.
#' @param angle `numeric` value indicating the rotation of ellipse.
#' @param segment NULL or `numeric` vector of two values indicating the
#'    start and end angles for the ellipse, prior to rotation.
#' @param arc.only `logical` indicating whether to draw the ellipse
#'    arc without connecting to the center of the ellipse. Set
#'    `arc.only=FALSE` when segment does not include the full circle,
#'    to draw only the wedge.
#' @param nv `numeric` the number of vertices around the center to draw.
#' @param deg `logical` indicating whether input `angle` and `segment`
#'    values are in degrees, or `deg=FALSE` for radians.
#' @param border,col,lty,lwd arguments passed to `graphics::polygon()`.
#' @param draw `logical` indicating whether to draw the ellipse.
#' @param ... additional arguments are passed to `graphics::polygon()`
#'    when `draw=TRUE`.
#'
#' @family jam igraph functions
#'
#' @examples
#' par("mar"=c(2, 2, 2, 2));
#' plot(NULL,
#'    type="n",
#'    xlim=c(-5, 20),
#'    ylim=c(-5, 18),
#'    ylab="", xlab="", bty="L",
#'    asp=1);
#' xy <- drawEllipse(
#'    x=c(1, 11, 11, 11),
#'    y=c(1, 11, 11, 11),
#'    a=c(5, 5, 5*1.5, 5),
#'    b=c(2, 2, 2*1.5, 2),
#'    angle=c(20, -15, -15, -15),
#'    segment=c(0, 360, 0, 120, 120, 240, 240, 360),
#'    arc.only=c(TRUE, FALSE, FALSE, TRUE),
#'    col=jamba::alpha2col(c("red", "gold", "dodgerblue", "darkorchid"), alpha=0.5),
#'    border=c("red", "gold", "dodgerblue", "darkorchid"),
#'    lwd=1,
#'    nv=99)
#' points(x=c(1, 11), y=c(1, 11), pch=20, cex=2)
#' jamba::drawLabels(x=c(12, 3, 13, 5),
#'    y=c(14, 10, 9, 2),
#'    labelCex=0.7,
#'    drawBox=FALSE,
#'    adjPreset=c("topright", "left", "bottomright", "top"),
#'    txt=c("0-120 degrees,\nangle=-15,\narc.only=TRUE",
#'       "120-240 degrees,\nangle=-15,\narc.only=TRUE,\nlarger radius",
#'       "240-360 degrees,\nangle=-15,\narc.only=FALSE",
#'       "angle=20"))
#'
#' @export
drawEllipse <- function(x,
                        y,
                        a=1,
                        b=1,
                        angle=0,
                        segment=NULL,
                        arc.only=TRUE,
                        nv=100,
                        deg=TRUE,
                        border=NULL,
                        col=NA,
                        lty=1,
                        lwd=1,
                        draw=TRUE,
                        ...){
  ## Purpose is to draw an ellipse
  if (length(deg) > 0 && any(deg %in% TRUE)) {
    deg <- TRUE
  } else {
    deg <- FALSE
  }
  
  if (length(segment) == 0) {
    if (length(deg) > 0 && any(deg %in% TRUE)) {
      segment <- c(0, 360);
    } else {
      segment <- c(0, pi*2);
    }
  }
  if (length(nv) == 0) {
    nv <- 100;
  } else if (length(nv) > 1) {
    nv <- head(nv, 1)
  }
  
  ## Fix various vector lengths
  y <- rep(y, length.out=length(x));
  a <- rep(a, length.out=length(x));
  b <- rep(b, length.out=length(x));
  col <- rep(col, length.out=length(x));
  border <- rep(border, length.out=length(x));
  
  ## if input is in degrees
  if (deg) {
    angle <- angle * pi/180;
    segment <- segment * pi/180;
  }
  segment <- rep(segment,
                 length.out=length(x) * 2);
  segment_seq <- seq(from=1, to=length(segment), by=2);
  segment1 <- segment[segment_seq];
  segment2 <- segment[segment_seq + 1];
  
  angle <- rep(angle,
               length.out=length(segment1));
  if (length(arc.only) == 0) {
    arc.only <- TRUE;
  }
  arc.only <- rep(arc.only,
                  length.out=length(segment1));
  if (length(segment1) == 1) {
    z <- seq(from=segment[1],
             to=segment[2],
             length=nv + 1);
    if (!arc.only) {
      z <- c(NA, z, NA, NA);
    }
    z_idx <- rep(1, length(z));
    z_angle <- rep(angle, length.out=length(z));
    z_cumsum <- length(z);
    z_lengths <- length(z);
  } else {
    z_list <- lapply(seq_along(segment1), function(i){
      j <- c(seq(from=segment1[i],
                 to=segment2[i],
                 length.out=nv + 1), NA);
      if (!arc.only[i]) {
        j <- c(NA, j, NA);
      }
      j
    })
    # jamba::printDebug("z_list:");print(z_list)
    z <- unlist(z_list);
    z_lengths <- lengths(z_list);
    z_idx <- rep(seq_along(z_list), z_lengths);
    z_angle <- rep(angle, z_lengths);
    z_cumsum <- cumsum(z_lengths);
  }
  xx <- a[z_idx] * cos(z);
  yy <- b[z_idx] * sin(z);
  alpha <- xyAngle(xx,
                   yy,
                   directed=TRUE,
                   deg=FALSE);
  rad <- sqrt(xx^2 + yy^2)
  xp <- rad * cos(alpha - z_angle) + x[z_idx];
  yp <- rad * sin(alpha - z_angle) + y[z_idx];
  if (any(!arc.only)) {
    which_wedge <- which(!arc.only);
    # jamba::printDebug("which_wedge: ", which_wedge);
    wedge_x <- x[which_wedge];
    wedge_y <- y[which_wedge];
    wedge_idx1 <- z_cumsum[which_wedge] - z_lengths[which_wedge] + 1;
    wedge_idx2 <- z_cumsum[which_wedge] - 1;
    # jamba::printDebug("wedge_idx1: ", wedge_idx1);
    # jamba::printDebug("wedge_idx2: ", wedge_idx2);
    xp[wedge_idx1] <- wedge_x;
    xp[wedge_idx2] <- wedge_x;
    yp[wedge_idx1] <- wedge_y;
    yp[wedge_idx2] <- wedge_y;
  }
  if (draw) {
    polygon(xp,
            yp,
            border=border,
            col=col,
            lty=lty,
            lwd=lwd,
            ...);
  }
  return(invisible(list(
    x=xp,
    y=yp,
    z=z
  )));
  invisible(list(x=xp,
                 y=yp));
}

shape.ellipse.plot <- function(coords,
                               v=NULL,
                               params){
  vertex.color <- params("vertex", "color");
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v];
  }
  vertex.frame.color <- params("vertex", "frame.color");
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v];
  }
  vertex.frame.width <- params("vertex", "frame.width");
  if (length(vertex.frame.width) == 0) {
    vertex.frame.width <- 1;
  }
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v];
  }
  vertex.size <- 1/200 * params("vertex", "size");
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v];
  }
  vertex.ellipse.ratio <- params("vertex", "ellipse.ratio");
  if (length(vertex.ellipse.ratio) == 0) {
    vertex.ellipse.ratio <- 2;
  }
  if (length(vertex.ellipse.ratio) != 1 && !is.null(v)) {
    vertex.ellipse.ratio <- vertex.ellipse.ratio[v];
  }
  
  # angle of rotation of the ellipse, clockwise in degrees
  vertex.ellipse.angle <- params("vertex", "ellipse.angle");
  if (length(vertex.ellipse.angle) == 0) {
    vertex.ellipse.angle <- 0;
  }
  if (length(vertex.ellipse.angle) != 1 && !is.null(v)) {
    vertex.ellipse.angle <- vertex.ellipse.angle[v];
  }
  
  drawEllipse(x=coords[,1],
              y=coords[,2],
              a=vertex.size,
              b=vertex.size/vertex.ellipse.ratio,
              col=vertex.color,
              border=vertex.frame.color,
              lwd=vertex.frame.width,
              angle=vertex.ellipse.angle,
              draw=TRUE);
}
