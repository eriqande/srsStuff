#' Process the output of srs_covar for PCA, etc.
#'
#' Compute a lot of stuff and returns a list
#' @param A the output list from `srs_covar()`
#' @param npp the number of principal components to retain in a tibble for plotting
#' and facetting over the different PCs.
#' @export
process_srs_covar_output <- function(A, npp = 5) {

  # first for each individual i, get the number of sites at which
  # that individual had a read and also another individual, j, averaged
  # over all j.  This tells us on average how many markers went into
  # compute the covariance for each individual.
  m <- A$M
  diag(m) <- 0
  ave_sites <- tibble(
    vcf_name = A$sample_names,
    ave_sites = rowMeans(m) * ncol(m) / (ncol(m) - 1)
  )

  # do the whole eigendecomposition on the standard covariance matrix.
  eig <- eigen(A$Cov)
  colnames(eig$vectors) <- sprintf("PC-%02d", 1:ncol(eig$vectors))


  pca_tib <- tibble::as_tibble(eig$vectors[,1:npp]) %>%
    mutate(vcf_name = A$sample_names) %>%
    select(vcf_name, everything())

  pca_long <- pca_tib %>%
    tidyr::gather(., key = "PC", "val", -vcf_name)

  # then expand a grid of the possible comparisons (ordered)
  pca_pairs <- expand.grid(vcf_name = pca_tib$vcf_name,
                           PCx = sprintf("PC-%02d", 1:npp),
                           PCy = sprintf("PC-%02d", 1:npp),
                           stringsAsFactors = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(., pca_long, by = c("vcf_name", "PCx" = "PC")) %>%
    dplyr::rename(val_x = val) %>%
    dplyr::left_join(pca_long, by = c("vcf_name", "PCy" = "PC")) %>%
    dplyr::rename(val_y = val)  %>%
    left_join(ave_sites, by = "vcf_name")


  # now, we also want to put each pairwise covariance value, and its associated
  # number of sites used in the calculation, into a tibble
  CM <- tibble(
    vcf_name_1 = rep(A$sample_names, times = nrow(A$Cov)),
    vcf_name_2 = rep(A$sample_names, each = ncol(A$Cov)),
    covar = as.vector(A$Cov),
    num_sites = as.vector(A$M),
  )

  # and now we are going to the the same sort of thing, but on the covariance
  # matrix computed by the corrlations estimated a la Weir and Goudet
#   WaG <- (A$IBS - A$Mt_S) / (1 - A$Mt_S)
#   eig <- eigen(WaG)
#   colnames(eig$vectors) <- sprintf("PC-%02d", 1:ncol(eig$vectors))
#
#
#   pca_tib <- as_tibble(eig$vectors[,1:6]) %>%
#     mutate(vcf_name = A$sample_names) %>%
#     select(vcf_name, everything())
#
#   pca_long <- pca_tib %>%
#     tidyr::gather(., key = "PC", "val", -vcf_name)
#
#   pcp <- expand.grid(vcf_name = pca_tib$vcf_name,
#                      PCx = sprintf("PC-%02d", 1:6),
#                      PCy = sprintf("PC-%02d", 1:6),
#                      stringsAsFactors = FALSE) %>%
#     tibble::as_tibble() %>%
#     dplyr::left_join(., pca_long, by = c("vcf_name", "PCx" = "PC")) %>%
#     dplyr::rename(WaG_x = val) %>%
#     dplyr::left_join(pca_long, by = c("vcf_name", "PCy" = "PC")) %>%
#     dplyr::rename(WaG_y = val)
#
#   # then join those both and return them
#   left_join(pcp, pca_pairs,
#             by = c("vcf_name", "PCx", "PCy"))


  list(
    pca_pairs = pca_pairs,
    covars_and_num_sites = CM,
    eig = eig
  )
}
