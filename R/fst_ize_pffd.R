
#' function that fst-izes the output from pairwise fst from dumpfile
#'
#' blah, blah
#' @param pffd  The output list from [pairwise_fst_from_dumpfile()].
#' @param treelist The named list that tells us who is in which populations.
#' We need this to know which population each index value in pffd represents.
#' @export
fst_ize_pffd <- function(pffd, treelist) {
  # now we get the popZbars
  popZbars <- pffd$popZ %>%
    as_tibble() %>%
    mutate(zbar_pop = z_sum / L_used) %>%
    select(idx, zbar_pop)

  # and we get the pairZbars
  pairZbars <- pffd$pairZ %>%
    as_tibble() %>%
    mutate(zbar_pair = pairSum / pairL) %>%
    select(pop1, pop2, zbar_pair)

  # and finally we can get the pairwise pop-specific fst's
  pwFs <- pairZbars %>%
    left_join(popZbars %>% rename(zbar_pop1 = zbar_pop), by = c("pop1" = "idx")) %>%
    left_join(popZbars %>% rename(zbar_pop2 = zbar_pop), by = c("pop2" = "idx")) %>%
    mutate(
      fst1 = (zbar_pop1 - zbar_pair) / (1 - zbar_pair),
      fst2 = (zbar_pop2 - zbar_pair) / (1 - zbar_pair)
    ) %>%
    mutate(
      Pop1 = names(treelist)[1:10][pop1],
      Pop2 = names(treelist)[1:10][pop2],
      ave_fst = (fst1 + fst2) / 2
    ) %>%
    select(Pop1, Pop2, fst1, fst2, ave_fst, everything()) %>%
    arrange(Pop1, ave_fst)

  pwFs
}
