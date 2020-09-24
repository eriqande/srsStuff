
## usethis namespace: start
#' @importFrom dplyr bind_rows case_when filter group_by_at left_join mutate pull rename select tally ungroup vars
#' @importFrom ggplot2 aes scale_fill_manual theme_void
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom magrittr %>%
#' @importFrom purrr map2 map_dbl
#' @importFrom Rcpp sourceCpp
#' @importFrom stats setNames
#' @importFrom tibble as_tibble tibble
#' @importFrom tidygraph as_tbl_graph graph_is_connected graph_is_tree node_is_leaf with_graph
#' @importFrom utils stack
## usethis namespace: end
NULL

## usethis namespace: start
#' @useDynLib srsStuff, .registration = TRUE
## usethis namespace: end
NULL

# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "daughter",
      "daughter_idx",
      "fst_boot_reps",
      "ind",
      "is_leaf",
      "label",
      "name",
      "node_type",
      "parent",
      "parent_idx",
      "text",
      "values"
    )
  )
}




