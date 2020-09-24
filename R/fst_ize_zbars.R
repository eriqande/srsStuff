#' Internal function to convert the z_bars to pop-specific Fst's
#'
#' Blah, blah
#' @param SI the list that gets output by srs_identity.
#' @param FET the fst_edges_tibble.  A tibble with columns parent daughter parent_idx and daughter_idx.
#' each daughter appears once.  The idxs must be the indexes of the columns in the
#' zbar output of SI that correspond to each parent or daughter node.
#' @keywords internal
fst_ize_zbars <- function(SI, FET) {
  FET %>%
    mutate(
      fst_boot_reps = map2(
        .x = parent_idx,
        .y = daughter_idx,
        .f = ~ (SI$zbar[, .y] - SI$zbar[, .x]) / (1 - SI$zbar[, .x])
      )
    ) %>%
    mutate(
      fst_est = map_dbl(fst_boot_reps, ~ .x[1]),
      min_boot = map_dbl(fst_boot_reps, ~ min(.x)),
      max_boot = map_dbl(fst_boot_reps, ~ max(.x))
    )
}
