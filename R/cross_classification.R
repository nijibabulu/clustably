#' @importFrom dplyr "%>%"
NULL


#' Align the clusters of two optionally downsampled cell clusterings
#'
#' This takes a greedy approach to aligning the identities of two
#' clustering partitions. Cells present in both sets are cross-
#' classified by their own labels. The pairs of labels with the
#' highest frequencies in descending order are assumed to be identical.
#' All cells in labels which are not present in the target are relabeled
#' with the default labeling.
#'
#' @param seuratQuery a seurat object to relabel
#' @param seuratTarget a seurat object with the target labels to find
#' @param defaultLabel label to give all classes not in the target
#'
#' @importFrom Seurat Idents
#'
#' @return The seuratQuery object with altered Idents
#' @export
align2Idents.Seurat <- function(seuratQuery, seuratTarget, defaultLabel="X") {
  newQueryIdents <- align2Idents(Idents(seuratQuery), Idents(seuratTarget), defaultLabel)
  seuratQuery[["aligned"]] <- newQueryIdents
  Idents(seuratQuery) <- "aligned"
  seuratQuery
}

#' Align the clusters of two downsampled cell labels
#'
#' This takes a greedy approach to aligning the identities of two
#' clustering partitions. Cells present in both sets are cross-
#' classified by their own labels. The pairs of labels with the
#' highest frequencies in descending order are assumed to be identical.
#' All cells in labels which are not present in the target are relabeled
#' with the default labeling.
#'
#' @param seuratQuery a seurat object to relabel
#' @param seuratTarget a seurat object with the target labels to find
#' @param defaultLabel label to give all classes not in the target
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange distinct mutate
#' @importFrom forcats fct_recode
#' @importFrom purrr set_names
#'
#' @return The seuratQuery object with altered Idents
#' @export
align2Idents <- function(query, target, defaultLabel="X") {
  aligned <-
    identsTable(target, query, names=c("target","query")) %>%
    as_tibble() %>%
    arrange(desc(n)) %>%
    distinct(query, .keep_all=T) %>%
    distinct(target, .keep_all=T)

  # rename the absent identities in the query to the default label
  absentIdents <- setdiff(query, aligned$query)
  if(length(absentIdents)) {
    absentIdentsMap <- set_names(absentIdents, length(absentIdents) %>% replicate(defaultLabel))
    query <- fct_recode(query, !!!absentIdentsMap)
  }

  targetIdentsMap <- set_names(aligned$query, aligned$target)
  query <- fct_recode(query, !!!targetIdentsMap)

  query
}


.mode <- function(xs) {
  uxs <- unique(na.omit(xs))
  uxs[which.max(tabulate(match(xs, uxs)))]
}

.confidence <- function(xs) {
  length(which(xs == .mode(xs))) / length(na.omit(xs))
}

#' Do the labeling of samples given a labeling function
#'
#' @param obj a seurat object
#' @param samples a set of cells to label
#' @param labelingFunc a function to label the cell subsets. should return a factor of labels named by the cells
#'                     This can be accomplished by just calling Idents() after clustering.
#' @param ncores a number of cores to use if multiprocessing is available. Uses the furrr library if greater than 1
#' @param verbose display a progress bar
#'
#' @return a list of labeling partitions on the samples
labelSamples <- function(obj, samples, labelingFunc, ncores, verbose) {
  if(ncores == 1) {
    pb <- progress_estimated(length(samples))
    .dolabel <- function(x) {
      if(verbose) { pb$tick()$print() }
      labelingFunc(subset(obj, cells=x))
    }
    labels <- map(samples, .dolabel)
    pb$stop()
  }
  else {
    plan(multiprocess, workers=ncores)
    labels <- future_map(samples,
                         ~labelingFunc(subset(obj, cells=.x)),
                         .progress=verbose)
    plan(sequential)
  }
  labels
}

#' add the cluster labels from an aligned tibble
#'
#' @param obj a seurat object
#' @param tbl a tibble containing columns with cell names and rows with labels
#' @param consensusName  store the consensus idents under this name
#' @param confidenceName store the confidence statistics under this name
#' @param freqName store the frequency of each labeling under this name, in the misc slot.
#' @param rawName store the raw table (the passed tbl) parameter under this name, in the misc slot.
#'
#' @importFrom purrr map_int map_chr map_dbl flatten_chr
#' @importFrom dplyr bind_rows
#' @importFrom tibble column_to_rownames add_column
#'
#' @return the Seurat object augmented with the consensus statistics
addClusterLabels <- function(obj, tbl, consensusName, confidenceName, freqName, rawName) {
  .mode <- function(xs) {
    uxs <- unique(na.omit(xs))
    uxs[which.max(tabulate(match(xs, uxs)))]
  }

  .confidence <- function(xs) {
    length(which(xs == .mode(xs))) / length(na.omit(xs))
  }

  .freq <- function(xs, labels) {
    map_int(labels, ~length(which(xs == .x)))
  }

  labels <- unique(tbl %>% flatten_chr() %>% na.omit())
  obj@misc[[freqName]] <- map(tbl, .freq, labels) %>% bind_rows() %>% add_column(rowname=labels) %>% column_to_rownames()
  obj@misc[[rawName]] <- tbl

  labelModes <- map_chr(tbl, .mode)
  obj[[consensusName]] <- factor(labelModes, levels=labelSort(unique(labelModes)))

  obj[[confidenceName]] <- map_dbl(tbl, .confidence)

  obj
}


# Get a bootstrap of the cluster labelings
#
# @param obj a seurat object. A subset of this object will be passed as the first parameter to the labelingFunc
# @param labelingFunc a function to label the cell subsets. should return a factor of labels named by the cells
#                     This can be accomplished by just calling Idents() after clustering.
# @param n number of times to apply bootstrap
# @param frac downsample fraction for bootstrap
# @param seed a random seed
# @param ncores do multicore processing if available
# @param verbose show progress
# @param consensusName store the consensus idents under this name
# @param bootstrapName store the consensus statistics under this name
# @param freqName store the frequency of each labeling under this name
# @param ... additional parameters for labelingFunc
#
# @return an updated seurat object
#
# @importFrom purrr map
#
# @export
bootstrapClustering <- function(obj, labelingFunc, n=10,  frac=.75, seed=42, ncores=1, verbose=FALSE,
                                consensusName="bootstrapConsensus",
                                bootstrapName="bootstrapConfidence",
                                freqName="bootstrapFrequency",
                                rawName="bootstrapResults",
                                ...) {
  set.seed(seed)

  cellSamples <- map(1:n, ~sample(Cells(obj), length(Cells(obj)), replace=TRUE))
  labels <- labelSamples(obj, cellSamples, labelingFunc, ncores, verbose)

  alignedTbl <- alignIdents(labels, Cells(obj))

  obj <- addClusterLabels(obj, alignedTbl, consensusName, bootstrapName, freqName, rawName)

  obj
}

#' Get a jackknife of the cluster labelings
#'
#' @param obj a seurat object. A subset of this object will be passed as the first parameter to the labelingFunc
#' @param labelingFunc a function to label the cell subsets. should return a factor of labels named by the cells
#'                     This can be accomplished by just calling Idents() after clustering.
#' @param n number of times to apply jackknife
#' @param frac downsample fraction for jackknife
#' @param seed a random seed
#' @param ncores do multicore processing if available
#' @param verbose show progress
#' @param consensusName store the consensus idents under this name
#' @param jackknifeName store the consensus statistics under this name
#' @param freqName store the frequency of each labeling under this name
#' @param ... additional parameters for labelingFunc
#'
#' @return an updated seurat object
#'
#' @importFrom purrr map
#'
#' @export
jackknifeClustering <- function(obj, labelingFunc, n=10,  frac=.75, seed=42, ncores=1, verbose=FALSE,
                                consensusName="jackknifeConsensus",
                                jackknifeName="jackknifeConfidence",
                                freqName="jackknifeFrequency",
                                rawName="jackknifeResults",
                                ...) {
  set.seed(seed)

  sampleSize = length(Cells(obj))*frac

  cellSamples <- map(1:n, ~sample(Cells(obj), sampleSize))
  labels <- labelSamples(obj, cellSamples, labelingFunc, ncores, verbose)

  encoding <- alignLabels(labels, verbose=verbose)
  recodedLabels <- recodeLabels(labels, encoding, Cells(obj))

  obj <- addClusterLabels(obj, recodedLabels, consensusName, jackknifeName, freqName, rawName)

  obj
}
