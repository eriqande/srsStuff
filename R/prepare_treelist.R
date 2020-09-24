
#' Prepare a tree of individuals/populations for population-specific Fst calculations
#'
#' This takes a list-based specification of a tree and it: checks it to make sure
#' it is correct, verifies that it is a tree, makes some plots of it so the
#' user can ensure they have input the data correctly, identifies the populations
#' (which are, by definition, one node level above the leaves), and then translates
#' all the portions of the tree into an integer-based specification that is
#' appropriate for the [srs_identity()] function.
#' @param treelist the desired tree, specified as a named list.  Each component
#' of this list is named for an internal node in the tree, and the contents of that
#' component is a character vector of the nodes (or leaves/individuals) that are the
#' daughters of that node. This tree is thought of as a rooted tree, with a topmost node
#' that is typically called "root", but could be named anything.  The first components
#' of the list must all the nodes representing populations, whose daughters are
#' exclusively individuals.  After that, the remaining nodes must be the other
#' internal nodes in the tree.
#' @param vcf_names  the names of the samples in the vcf file in the order that their
#' columns appear in the vcf file.  These don't have to be exactly the names as they
#' appear in the VCF file.  The important part is that they correspond to the individual
#' names in treelist and that they are ordered as the individuals appear as columns
#' in the VCF file. Note that in this current implementation, the columns of the VCF
#' file must not include any individuals other than those that are in treelist. I might
#' try to change that later, but have no guarantees.
#' @export
prepare_treelist <- function(treelist, vcf_names = NULL) {

  ret <- list()

  # check to make sure none of the nodes appear twice as a daughter
  dupie_daughters <- unlist(treelist)[duplicated(unlist(treelist))]
  if(length(dupie_daughters) > 0) {
    stop(
      "Problem, treelist has these duplicate daughter entries: ",
      paste(dupie_daughters, collapse = ", ")
    )
  }

  # check to make sure that none of the internal nodes are given twice
  dupie_internal <- names(treelist)[duplicated(names(treelist))]
  if(length(dupie_internal) > 0) {
    stop(
      "Problem, treelist has these duplicate internal node entries: ",
      paste(dupie_internal, collapse = ", ")
    )
  }


  # first, turn the list into a tibble of edges defined by the nodes at
  # which they start and stop
  # make a data frame of edge start and stops
  d <- lapply(names(treelist), function(n) {
    tibble(parent = n, daughters = treelist[[n]])
  }) %>%
    bind_rows()

  # make a tidygraph from that
  tgraph <- as_tbl_graph(d)

  # verify it is connected
  if (!with_graph(tgraph, graph_is_connected())) {
    stop("Sorry, the elements of treelist are not connected in a single tree.")
  }

  # verify it is a tree
  if (!with_graph(tgraph, graph_is_tree())) {
    stop("Sorry, treelist is not a tree.")
  }





  # Now, identify the "population" nodes.  These are the ones
  # in which all the daughters are leaves.  At some point I might want
  # to consider a situation where the daughters of a node are a mix of
  # leaves and internal nodes.  But the C++ code is not set up for that
  # at the moment. So, I enforce this defn of a population node currently.
  tg2 <- tgraph %>%
    mutate(is_leaf = node_is_leaf())

  leaves <- tg2 %>%
    filter(is_leaf) %>%
    pull(name)

  # first, throw an error if any node has a mix of leaves and internal nodes
  # under it
  mixed_nodes_lgl <- lapply(
    treelist,
    function(x) {
      any(x %in% leaves) && any(!(x %in% leaves))
    }
  ) %>%
    unlist()
  if (sum(mixed_nodes_lgl) > 0) {
    stop(
      "Sorry. In the current implementation each node's daughters must all ",
      "be leaves (individuals) or internal nodes.  The following nodes are mixed: ",
      paste(names(mixed_nodes_lgl)[mixed_nodes_lgl], collapse = ", ")
    )
  }


  # here is a logical vector giving the populations
  pop_nodes_lgl <- lapply(
    treelist,
    function(x) {
      all(x %in% leaves)
    }
  ) %>%
    unlist()

  # check to make sure that these population nodes all come first
  if (!(all(pop_nodes_lgl[1:sum(pop_nodes_lgl)] == TRUE))) {
    stop(
      "Found ",
      sum(pop_nodes_lgl),
      " population nodes, but the first ",
      sum(pop_nodes_lgl),
      " in treelist are not all population nodes."
      )
  }

  #### At this juncture, let's make some plots of this tree ####
  tg3 <- tg2 %>%
    mutate(
      label = ifelse(is_leaf, "", name),
      node_type = case_when(
        is_leaf == TRUE ~ "individual",
        name %in% names(treelist)[pop_nodes_lgl] ~ "population",
        TRUE ~ "ancestral/internal"
      )
    )

  # a quick function to add nodes and edges to the ggraph
  add_nodes_etc <- function(g) {
    g +
      geom_edge_link() +
      geom_node_point(aes(fill = node_type), shape = 21, stroke = 0.3) +
      scale_fill_manual(
        values = c(
          `ancestral/internal` = "blue",
          individual = "gray",
          population = "orange"
        )
      ) +
      theme_void()
  }


  # plot the tree
  g1 <- add_nodes_etc(ggraph(tg3, layout = "fr"))
  # make a dendrogram version too
  d1 <- add_nodes_etc(ggraph(tg3, layout = "dendrogram"))


  # this function adds some labels to the non-population internal nodes
  add_node_labels <- function(g) {
    g +
      geom_node_text(aes(label = label),   colour = "yellow") +
      geom_node_text(aes(label = label),  fontface = "bold", colour = "black")
  }

  plots <- list()
  plots$tree_without_names <- g1
  plots$tree_with_names <- add_node_labels(g1)
  plots$dendrogram_without_names <- d1
  plots$dendrogram_with_names <- add_node_labels(d1)
  ret$plots <- plots

  #### NOW, PREPARE INPUT FOR srs_identity() ####
  if (!is.null(vcf_names)) {

    # first, check to make sure that each individual in treelist appears in there
    indivs <- tg3 %>%
      filter(is_leaf) %>%
      pull(name)

    missers <- setdiff(indivs, vcf_names)
    if (length(missers) > 0) {
      stop(
        "The following individuals in treelist do not appear in vcf_names: ",
        paste(missers, collapse = ", ")
      )
    }
    missers <- setdiff(vcf_names, indivs)
    if (length(missers) > 0) {
      stop(
        "The following individuals in vcf_names do not appear in vcf_names: ",
        paste(missers, collapse = ", "),
        "  Big apologies, but every individual in the extracted VCF must be represented in treelist"
      )
    }

    # record the number of populations and get their indexes
    num_pops <- sum(pop_nodes_lgl)

    # record the number of internal nodes
    num_internal_nodes <- length(treelist)


    # get the 0-based index of the population each individual is in,
    # **in the order of individuals in the VCF file**
    internal_0_indexes <- 0:(num_internal_nodes - 1)
    names(internal_0_indexes) <- names(treelist)

    tmp <- as_tibble(stack(treelist)) %>%
      mutate(ind = as.character(ind)) %>%
      mutate(idx = internal_0_indexes[ind])
    pop_0_indexes <- tmp$idx
    names(pop_0_indexes) <- tmp$values

    pops_of_indivs <- unname(pop_0_indexes[vcf_names])

    # finally, for the non-population internal nodes, let's make the list
    # of daughters for them all.
    internal_node_daughters <- lapply(
      treelist[(num_pops + 1):length(treelist)],
      function(x) unname(internal_0_indexes[x])
    )

    # record all that into a list that we can output
    srs_id_input <- list()
    srs_id_input$pops_of_indivs <- pops_of_indivs
    srs_id_input$num_pops <- num_pops
    srs_id_input$num_internal_nodes <- num_internal_nodes
    srs_id_input$non_pop_internal_node_daughters <- internal_node_daughters

    #### Now, prepare some information to be used to extract the results from srs_identity() into something meaningful ####
    # we need a listing of the edges keyed out by daughter index (base-1 this time).
    fst_edges_tibble <- as_tibble(stack(treelist[(num_pops + 1):length(treelist)])) %>%
      mutate(ind = as.character(ind)) %>%
      rename(
        parent = ind,
        daughter = values
      ) %>%
      select(
        parent,
        daughter
      ) %>%
      mutate(
        parent_idx = internal_0_indexes[parent] + 1,
        daughter_idx = internal_0_indexes[daughter] + 1,
      )

    post_proc_stuff <- list()
    post_proc_stuff$fst_edges_tibble <- fst_edges_tibble

    ret$srs_id_input <- srs_id_input
    ret$srs_post_proc_input <- post_proc_stuff

  } # closes: if (!is.null(vcf_names))


  return(ret)
}
