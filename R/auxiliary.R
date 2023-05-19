fill.na <- function(x, i = 5) {
  if (is.na(x)[i]) {
    return(round(mean(x, na.rm = TRUE), 0))
  }
  else {
    return(round(x[i], 0))
  }
}

draw.circle <- function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1,
                         density = NULL, angle = 45, lwd = 1) {
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- plotrix::getYmult()
  angle.inc <- 2 * pi/nv
  angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
  if (length(col) < length(radius))
    col <- rep(col, length.out = length(radius))
  for (circle in 1:length(radius)) {
    xv <- cos(angles) * radius[circle] + x
    # yv <- sin(angles) * radius[circle] * ymult + y
    yv <- sin(angles) * radius[circle] + y
    polygon(xv, yv, border = border, col = col[circle], lty = lty,
            density = density, angle = angle, lwd = lwd)
  }
  invisible(list(x = xv, y = yv))
}

slotApply <- function(x,FUN,...){
  cl <- class(x)
  result <- list()
  for(i in slotNames(cl)){
    result[[i]] <- FUN(slot(x,i),...)
  }
  result
}

slotToList <- function(x){
  returnlist <- list()
  namesslot <- slotNames(x)
  for(cur_slot in namesslot)
    returnlist[[cur_slot]] <- slot(x, name = cur_slot)
  returnlist
}

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
