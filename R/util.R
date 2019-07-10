#' labelSort
#
#' Sort labels numerically if possible and then alphabetically.
#
#' @param xs a character vector
#
#' @return a sorted character vector
#'
#' @export
labelSort <- function(xs) {
  is.num <- function(y) !is.na(suppressWarnings(as.numeric(y)))
  as.chrnum <- function(y) ifelse(is.num(y), as.numeric(y), y)
  '>.ic' <<- function(a,b) {
    if(is.num(a) == is.num(b))  return(as.chrnum(a) > as.chrnum(b))
    else return(is.num(a))
  }
  '==.ic' <<- function(a, b) ifelse(a > b || b > a, FALSE, TRUE)
  '<.ic' <<-  function(a, b) b > a

  '[.ic' <<- function(x, i) {
    class(x) <- "character"
    x <- x[i]
    class(x) <- "ic"
    x
  }

  class(xs) <- "ic"
  sort(xs)
}

#' Retrieve a table of cross-classification table of cell identities
#'
#' @param x a named factor
#' @param y a named factor
#'
#' @return an R table
#' @export
identsTable <- function(x, y, names=c("x","y")) {
  commonCells <- intersect(names(x), names(y))
  table(x[commonCells],  y[commonCells], dnn=names)
}
