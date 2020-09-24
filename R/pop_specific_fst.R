

#' Compute population-specific Fst's for populations arranged on a hierarchical tree
#'
#' Need to write more
#' @param in_list a list of inputs like that returned from `prepare_treelist()`.
#' @param AD_file a file with the allele depths at the markers. The columns in this
#' file must be ordered according to the vcf_names passed to `prepare_treelist()`.
#' @param boot_reps the number of Poisson bootstrap samples to compute Fst for.
#' @export
pop_specific_fst <- function(
  in_list,
  AD_file,
  boot_reps = 100
  ) {

  L <- in_list$srs_id_input
  boot_reps <- as.integer(boot_reps)
  AD_file_N <- normalizePath(AD_file)

  SID <- srs_identity(
    file = AD_file_N,
    pops_of_indivs = as.integer(L$pops_of_indivs),
    num_pops = as.integer(L$num_pops),
    num_internal_nodes = as.integer(L$num_internal_nodes),
    daughters = L$non_pop_internal_node_daughters,
    BootReps = boot_reps
    )

  Fsts <- fst_ize_zbars(
    SI = SID,
    FET = in_list$srs_post_proc_input$fst_edges_tibble
  )

  Fsts

}
