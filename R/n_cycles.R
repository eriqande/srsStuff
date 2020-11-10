#' Find directed cycles of a given length in a directed graph
#'
#' This is intended to be used iteratively: first find all the cycles of
#' length 2.  Then remove any edges that involve any of the nodes in those
#' two-cycles.  Then see if there are any three-cycles, etc.
#' @param G the directed graph passed in as a list with names as the parents
#' and daughters as the values.  The graph must be such that each
#' parent has only a single daughter.  Though a node can be the daughter
#' of multiple parents.
#' @param n the length of the cycles that should be found
#' @return Returns a list with the following elements:
#' - `length` holds the values of `n`
#' - `cycles` is a list of the cycles in which each cycle is
#' just the names of the nodes in cycle order, starting with the first one,
#' lexicographically.
#' - `remaining` is a list of edges remaining after all of the nodes involved
#' in the currently-found cycles have been removed.
#' @export
n_cycles <- function(G, n) {

  # if requesting two-cycles, that is pretty quick and we can do it like this
  if(n == 2) {
    two_cycles <- lapply(names(pp), function(n) {
      if(pp[[pp[[n]]]] == n)
        return(sort(c(n, pp[[n]])))
      else
        return(character(0))
    }) %>%
      unique()

    # then get rid of the empty ones:
    two_cycles <- two_cycles[!(sapply(two_cycles, length) == 0)]

    # then remove the nodes in the two cycles (both as parents and daughters)
    tcf <- flatten(two_cycles)
    rpp <- G
    rpp <- rpp[!(names(rpp) %in% tcf)]  # toss out parents in the two-cycles
    rpp <- rpp[!sapply(rpp, function(x) x %in% tcf)]  # toss out daughters in the two-cycles

    # finally, return the list
    ret <- list(
      length = n,
      cycles = two_cycles,
      remaining = rpp
    )
    return(ret)
  }

  # down here we deal with the longer cycles
}
