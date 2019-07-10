
#' Perform a global alignment of labelings on a set of objects
#'
#' Greedily group cell labelings and recode them into a single tibble
#'
#' @param labels a list of factors with labelings
#' @param cells a character vector of all cells ids in the experiment
#'
#' @importFrom tibble as_tibble enframe
#' @importFrom dplyr any_vars arrange bind_rows desc filter_all mutate mutate_if right_join select
#' @importFrom purrr map2 map_dfc set_names
#' @importFrom forcats fct_recode
#' @importFrom tidyr crossing spread
#'
#' @return a tibble containing the recoded labelings
alignIdents <- function(labels, cells) {
  n <- length(labels)
  labelTab <- table(labels, dnn=1:n) %>% as_tibble() %>% arrange(desc(n))
  encodingTab <- labelTab %>% select(-n) %>% map_dfc(~ifelse(duplicated(.x), NA, .x)) %>% filter_all(any_vars(!is.na(.)))
  recodedLabels <- map2(labels, encodingTab, ~fct_recode(.x, !!!set_names(.y, 1:length(.y)) %>% na.omit()))

  alignedTbl <-
    seq(n) %>%
    map(~enframe(recodedLabels[[.x]], name="cell", value="ident") %>%
          mutate_if(is.factor, as.character) %>%
          mutate(n=.x)) %>%
    bind_rows()  %>%
    right_join(crossing(n=seq(n), cell=cells)) %>%
    spread(cell,ident) %>%
    select(-n)

  alignedTbl
}

#' Recode labels by an encoding of the labels after alignment
#'
#' @param labels the original labels given to to the cells
#' @param encoding the encoding given by alignLabels
#' @param cells all the cell names in the expeiriment
#'
#' @importFrom purrr map map2
#' @importFrom dplyr bind_rows right_join select
#' @importFrom tidyr spread
#' @importFrom forcats fct_recode
#'
#' @return a tbl with columns of cell names and rows of each of the replicates
#' @export
recodeLabels <- function(labels, encoding, cells) {
  n <- length(labels)
  recodedLabels <- map2(labels,
                        encoding[as.character(1:length(encoding))],
                        ~fct_recode(.x, !!!set_names(.y, 1:length(.y)) %>% na.omit()))

  recodedTbl <-
    seq(n) %>%
    map(~enframe(recodedLabels[[.x]], name="cell", value="ident") %>%
          mutate_if(is.factor, as.character) %>%
          mutate(n=.x)) %>%
    bind_rows()  %>%
    right_join(crossing(n=seq(n), cell=cells), by=c("cell","n")) %>%
    spread(cell,ident) %>%
    select(-n)

  recodedTbl
}


# Find the top path among the pairs. Walk through each of the top pairs
# and add it if neither replicate has been seen or the one that has been
# seen has the same cluster label (i.e. do not add an incompatible link)
# this is roughly equivalent to a greedy max path through the matrix.
#
# @param info a list containing pairs - a tibble of inter-replicate label pairs sorted by frequencies
#                           and encoding - the current encoding of the replicates and pairs
# @para n the number of replicates in the data set
#
#' @importFrom purrr map2
#' @importFrom dplyr filter inner_join anti_join full_join group_by_all select summarize bind_rows
#' @importFrom tibble tibble
#' @importFrom tidyr separate
topPath <- function(info, n) {

  pairs <- info$pairs
  encoding = tibble(rep=c(pairs[1,]$x.rep, pairs[1,]$y.rep),
                    clust=c(pairs[1,]$x.clust, pairs[1,]$y.clust))
  for(i in 2:nrow(pairs)) {
    with(pairs[i,], {
          seen.x = encoding %>% filter(rep == x.rep)
          seen.y = encoding %>% filter(rep == y.rep)
          if(!(nrow(seen.x) && nrow(seen.y))  &&
             !(nrow(seen.x) && seen.x$clust != x.clust) &&
             !(nrow(seen.y) && seen.y$clust != y.clust)) {
              encoding <<-
                bind_rows(encoding,
                          tibble(rep=c(x.rep,y.rep),
                                 clust=c(x.clust,y.clust))) %>%
                distinct()
             }
        })
    if(nrow(encoding) == n) { break }
  }
  encodingRow <- encoding %>% full_join(tibble(rep=as.character(1:n)), by="rep") %>% spread(rep, clust)
  info$encoding <- bind_rows(info$encoding,  encodingRow)
  info$pairs <- info$pairs %>%
    anti_join(encoding, by=c("x.rep"="rep", "x.clust"="clust")) %>%
    anti_join(encoding, by=c("y.rep"="rep", "y.clust"="clust"))
  info
}


#' Perform a gloabl alignment of replicate labels.
#'
#' Greedy search for the most likely cluster relationship between replicates.
#'
#' @param labels a list of character vectors of labels named by their cell ID
#' @param verbose show progress
#'
#' @return a tibble of encodings of the clusters
#'
#' @importFrom purrr map2
#' @importFrom dplyr inner_join filter summarize progress_estimated
#' @importFrom tidyr separate
#' @importFrom stringr str_c
#'
#' @export
alignLabels <- function(labels, verbose=FALSE) {
  n <- length(labels)

  # recode the labels into a 2-column tbl with the cellID and the
  # labeling encoded as replicate_value
  labelsTbl <- map2(1:n, labels,
                    ~str_c(as.character(.x), .y, sep="_") %>%
                      set_names(names(.y)) %>%
                      enframe())  %>%
                bind_rows()

  # create a cross-classification tibble which is sorted by the number
  # of overlapping annotations between two cells and sort them in
  # descending order of
  pairs <- inner_join(labelsTbl, labelsTbl, by="name") %>%
    filter(value.x != value.y) %>% select(-name) %>% group_by_all() %>%
    summarize(n=n()) %>% arrange(desc(n)) %>% select(-n) %>%
    separate(value.x, c("x.rep", "x.clust"), "_") %>%
    separate(value.y, c("y.rep", "y.clust"), "_") %>%
    filter(x.rep != y.rep)

  if(verbose) { pb <- progress_estimated(nrow(pairs)) }
  info <- list(pairs=pairs, encoding=tibble())
  repeat {
    info <- topPath(info, n)
    if(verbose) {
      pb$i <- pb$n-nrow(info$pairs)
      pb$print()
    }
    if(nrow(info$pairs) == 0) { break }
  }
  if(verbose) { pb$stop() }

  info$encoding
}
