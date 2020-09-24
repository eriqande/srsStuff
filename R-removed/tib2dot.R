#' given a tibble specification of the pop-hierarchy, make a dot file of the tree
#'
#' @param tib the tibble
#' @param outfile the file to write it out.  By default it is tree.dot
#' @param num_gene_copies for each individual at the bottom of the tree,
#' how many gene copies should come off of him/her
#' @param int_fontsize size of font for labels on internal nodes.  User
#' must fiddle with this to get it right.
#' @param indiv_size size of the little circles at the bottom that denote individuals.
#' User must fiddle with this to get it right.
#' @export
tib2dot <- function(tib,
                    outfile = "tree.dot",
                    names_as_numbers = FALSE,
                    names_on_indivs = FALSE,
                    num_gene_copies = 2,
                    int_fontsize = 50,
                    indiv_size = 0.5,
                    edge_extras = NULL
                    ) {

  # start off writing the preamble
  cat("digraph G{\n\tratio=fill;\n\tnodesep=0.15;\n\tsize=\" 10, 7.5 \";\n\torientation=landscape;\n", file = outfile)

  # now, we write out the edge specs for all the edges
  cols <- ncol(tib)
  for(i in cols:2) {
    keepies <- tib[[i]] != tib[[i-1]]
    tmp <- tib %>%
      filter(keepies) %>%
      group_by_at(vars(i, i-1)) %>%
      tally() %>%
      ungroup() %>%
      setNames(c("parent", "daughter", "n")) %>%
      select(parent, daughter)

    if(!is.null(edge_extras)) {
      tmp2 <- tmp %>%
        left_join(edge_extras, by = c("parent", "daughter")) %>%
        mutate(text = ifelse(is.na(text), "", text))

      sprintf("\t\"%s\" -> \"%s\" [dir=none%s];\n", as.character(tmp2[[1]]), as.character(tmp2[[2]]), tmp2$text) %>%
        cat(file = outfile, append = TRUE)
    } else {
      sprintf("\t\"%s\" -> \"%s\" [dir=none];\n", as.character(tmp[[1]]), as.character(tmp[[2]])) %>%
        cat(file = outfile, append = TRUE)
    }
  }

  # and now make the internal nodes big, by increasing the font size
  setdiff(levels(tib[[1]]), tib[[1]]) %>%
    sprintf("\t\"%s\" [fontsize=%f, label=\"%s\"];\n", ., int_fontsize, .) %>%
    cat(file = outfile, append = TRUE)


  # now we write out a cluster that includes the individuals and the gene copies
  # which will put them all on one level and let us manipulate their ranks

  cat("\nsubgraph clusterpops {\n\tranksep=0.5;\n\tcolor=invis;\n", file = outfile, append = TRUE)


  # now, make the individuals small open circles
  sprintf("\t%s [shape=circle, height=%f, label=\"\"];\n", tib[[1]], indiv_size) %>%
    cat(file = outfile, append = TRUE)

  # and we add to that the edge-specs for the gene copies
  tmp <- list()
  for(i in 1:num_gene_copies) {
    tmp[[i]] <- paste(tib[[1]], i, sep = "_")
  }
  names(tmp) <- 1:num_gene_copies
  gtib <- as_tibble(tmp)

  # first do the nodes for the gene copies
  for(i in 1:ncol(gtib)) {
    sprintf("\t%s [shape=point, label=\"\"];\n", gtib[[i]]) %>%
      cat(file = outfile, append = TRUE)
  }
  # and then the edges leading to the gene copies
  for(i in 1:ncol(gtib)) {
    sprintf("\t\"%s\" -> \"%s\" [dir=none];\n", as.character(tib[[1]]), as.character(gtib[[i]])) %>%
      cat(file = outfile, append = TRUE)
  }

  cat("\n}", file = outfile, append = TRUE)  # add a closing brace for the cluster



  cat("\n}", file = outfile, append = TRUE)  # closing brace on the whole graph
}
