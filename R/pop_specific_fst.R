

#' Compute population-specific Fst's for populations arranged on a hierarchical tree
#'
#' Need to write more
#' @param in_list a list of inputs like that returned from `prepare_treelist()`.
#' @param AD_file a file with the allele depths at the markers. The columns in this
#' file must be ordered according to the vcf_names passed to `prepare_treelist()`.
#' @param boot_reps the number of Poisson bootstrap samples to compute Fst for.
#' @param dump_pop_Ys  if TRUE, then this will write out the allele counts at each
#' locus at each population.
#' @param dump_file  the path to the file to dump the Ys into.  By default, this
#' is set to a temporary file by R.  But, if you don't want it to disappear between
#' R sessions, you can put is somewhere else.
#' @export
pop_specific_fst <- function(
  in_list,
  AD_file,
  boot_reps = 100,
  dump_pop_Ys = FALSE,
  dump_file = tempfile()
  ) {

  dump_flag <- 0L
  if(dump_pop_Ys) {
    dump_flag <- 1L
    message("Dumping binary file of Ys to: ", dump_file)
  }

  L <- in_list$srs_id_input
  boot_reps <- as.integer(boot_reps)
  AD_file_N <- normalizePath(AD_file)

  SID <- srs_identity(
    file = AD_file_N,
    pops_of_indivs = as.integer(L$pops_of_indivs),
    num_pops = as.integer(L$num_pops),
    num_internal_nodes = as.integer(L$num_internal_nodes),
    daughters = L$non_pop_internal_node_daughters,
    BootReps = boot_reps,
    dump_pop_Ys = dump_flag,
    popYdump_file = dump_file
    )

  Fsts <- fst_ize_zbars(
    SI = SID,
    FET = in_list$srs_post_proc_input$fst_edges_tibble
  )

  list(
    Fsts = Fsts,
    srs_id_output = SID,
    in_list = in_list,
    dumpfile = dump_file
  )
}
