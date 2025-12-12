box::use(
  checkmate[assert_atomic_vector],
  chk = checkmate
)


#' Replace values of a vector based on a mapping scheme.
#'
#' @description
#' This is basically a copy or re-implementation of the \code{plyr::mapvalues()} function.
#' It takes an input vector \code{x} as well as two vectors \code{from} and
#' \code{to} that act as a dictionary for the mapping.
#'
#' @source
#' https://github.com/DavidJesse21/data.table.extras/blob/main/R/mapvalues.R
#'
#' @param x (`any`)\cr
#'   A vector.
#' @param from (`any`)\cr
#'   A vector specifying the values that should be mapped.
#' @param to (`any`)\cr
#'   A vector specifying the values to map to.
#'
#' @return (`any`)\cr
#'   A vector. It may be of a different type than the input vector `x`,
#'   depending on your inputs for `from` and `to`.
#'
#' @export
mapvalues = function(x, from, to) {
  assert_atomic_vector(x, min.len = 1L)
  assert_atomic_vector(from, min.len = 1L, unique = TRUE)
  assert_atomic_vector(to, min.len = 1L)
  
  if (!same_length(from, to)) {
    stop("Vectors `from` and `to` must have the same length.")
  }
  
  if (!all(unique(x) %in% from)) {
    stop("`from` must cover all unique values in `x`.")
  }
  
  
  if (is.factor(x)) {
    levels(x) = mapvalues(levels(x), from, to)
    return(x)
  }
  
  mapidx = match(x, from)
  return(to[mapidx])
}


same_length = function(x, y) {
  length(x) == length(y)
}



#' Check if `x` is a binary vector.
checkBinary = function(x) {
  check1 = chk$test_integerish(x)
  check2 = chk$test_set_equal(unique(x), 0:1)
  res = check1 && check2
  
  if (isFALSE(res)) {
    return("Must be binary (0 and 1)")
  } else {
    return(TRUE)
  }
}

#' @export
assert_binary = chk$makeAssertionFunction(checkBinary)


