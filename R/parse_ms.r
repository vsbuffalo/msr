
make_flag <- function(flag, vals, strict=FALSE) {
  # make a unix CL flag by appending one or two dashes. underscore turned into -
  dash <- "-"
  if (nchar(flag) > 1 && flag != 'seed')
    dash <- "--"
  # if a logical value is passed, e.g. T=TRUE, turn that into -T
  if (is.logical(vals)) {
    if (vals)
     return(sprintf("%s%s", dash, flag))
    # else, flag is off, don't include 
  }
  # if a underscore is present (e.g. for mspms's -random_seed) convert to dash
  flag <- gsub('_', '-', flag, fixed=TRUE) 
  sprintf("%s%s %s", dash, flag, paste(format(vals, scientific=FALSE), collapse=" "))
}

strict_make_flag <- function(flag, vals) {
  flag <- gsub('.', '-', flag, fixed=TRUE)
  sprintf("%s %s", flag, paste(format(vals, scientific=FALSE), collapse=" "))
}

make_args <- function(..., strict=FALSE) {
  args <- list(...)
  if (!strict)
    return(paste(Map(make_flag, names(args), args), collapse=" "))
  else
    return(paste(Map(strict_make_flag, names(args), args), collapse=" "))
}

#' Call MS from R
#'
#' Call Hudon's MS from R, returning a string of results (that can be parsed
#' with parse_ms()). Other coalescent simulators are supported as long as they
#' output results in MS format; specify with the argument \code{ms='mspms'}
#' (in this example, running msprime here). No argument checking occurs in this 
#' function; arguments passed to \code{...} are converted to command line 
#' arguments, e.g. \code{r=c(4, 1000)} is converted to \code{-r 4 1000} and 
#' \code{mutation.rate=1e-5} converted to \code{--mutation-rate 1e-5}.
#'
#' @param nsam number of samples (gametes) to draw
#' @param howmany how many replicates to run
#' @param cmd the command to pass to MS 
#' @param ... command line arguments as function arguments (if not using \code{cmd})
#' @param strict don't try to to guess flags by length of flag name; ".f 3" converted to -f 3, "..flag 3" converted to --flag 3
#' @param ms the ms-like executable to run, must be in \code{$PATH} or path to executable
#' 
#' @export
call_ms <- function(nsam, howmany, cmd=NULL, ..., ms="ms", strict=FALSE, verbose=TRUE) {
  func_args <- make_args(...)
  if (length(func_args) > 0 && !is.null(cmd))
    stop("specify string command line arguments through 'cmd' or arguments through '...', not both")
  if (length(func_args) > 0)
    cmd <- make_args(..., strict=strict)
  ms_cmd <- sprintf("%s %d %d %s", ms, nsam, howmany, cmd)
  if (verbose)
    message(sprintf("command: %s\n", ms_cmd))
  system(ms_cmd, intern=TRUE)
}

#' Call MS and Parse Output from R
#'
#' Call Hudon's MS from R, returning a tibble of results that has been parsed
#' by parse_ms(). This is simply a wrapper for call_ms() and parse_ms()
#'
#' @param nsam number of samples (gametes) to draw
#' @param howmany how many replicates to run
#' @param cmd the command to pass to MS 
#' @param ... command line arguments as function arguments (if not using \code{cmd})
#' @param ms the ms-like executable to run, must be in \code{$PATH} or path to executable
#
#' 
#' @export
ms <- function(nsam, howmany, cmd=NULL, ..., ms="ms", verbose=TRUE) {
  parse_ms(call_ms(nsam=nsam, howmany=howmany, cmd=cmd, ..., ms=ms, 
                   verbose=verbose))
}


#' A "tee" (in Unix sense) that writes MS's results to file
#'
#' @param x MS results @param con to save output to
#' @param pass if TRUE, pass the input directly to output for pipe operations
#' 
#' This writes a connection \code{con}, then return results (so can be used in
#' pipe).
#'
#' @export
#' @examples
#' # write MS results to 'test.ms', then pass down the pipeline
#' \dontrun{
#' res <- call_ms(10, 500, "-t 5") %>% write_tee('test.ms') %>% 
#'              parse_ms() %>% mutate(pi=map_dbl(gametes, pi))
#' }
write_tee <- function(x, con, pass=TRUE) {
  writeLines(x, con=con)
  if (pass)
    return(x)
}

#' Parse MS's key/value pairs, e.g. segsites and positions
#' returning a list of key/vals (where vals can be list too)
#' @keywords internal
parse_keyvals <- function(x) {
  keyvals <- gsub("(\\w+):\\s+(.*)", "\\1;;\\2", x, perl=TRUE)
  tmp <- strsplit(keyvals, ";;")[[1]]
  key <- tmp[1]
  vals <- as.numeric(strsplit(tmp[2], "\\s+")[[1]])
  if (length(vals) == 1)
    return(setNames(tibble::tibble(vals), key))
  else
    return(setNames(tibble::tibble(list(vals)), key))
}


#' Convert gaemtes' alleles intro matrix
#' @keywords internal
sites_matrix <- function(x) {
  do.call(rbind, lapply(x, function(y) as.integer(unlist(strsplit(y, "")))))
}

has_tree <- function(x) {
  # x are lines for a sim replicate. If there's a tree, second line
  # begins with (
  regexpr("^(\\(|\\[)", x[2]) != -1
}

has_many_trees <- function(x) regexpr("^\\[", x[2]) != -1

extract_recomb_trees <- function(x) {
  chunks <- strsplit(sub("\\[(\\d+)\\](.+)", "\\1;;;;\\2", x, perl=TRUE), ";;;;")
  lens <- as.integer(sapply(chunks, '[', 1))
  trees <- sapply(chunks, '[', 2)
  list(tibble::tibble(lens=lens, tree=trees))
}

#' Tidy a single simulation result from MS
#' @keywords internal
tidy_sim <- function(x) {
  # first element is the delimiter "//", last element is blank line (except for
  # last sim) 
  stopifnot(x[1] == "//")
  i <- 2
  keyval_i <- grep("^[a-z]+", x)
  if (has_tree(x)) {
    # for sims with recombination, we have multiple trees until segsites
    if (has_many_trees(x)) {
      # parse each tree, extracting the segment length (bp) and tree
      tree <- extract_recomb_trees(x[i:(min(keyval_i)-1)])
    } else {
      tree <- x[i]  # just one tree per locus, e.g. no recomb
    }
  }

  keyval_lines <- x[keyval_i]

  # parse keyval lines, and start creating the master tibble
  out <- dplyr::bind_cols(lapply(keyval_lines, parse_keyvals))

  stopifnot(nrow(out) == 1)
  if (out$segsites[1] > 0) {
    # extract gamete lines, remove empty line if there, convert gametes to matrices
    gametes <- x[(max(keyval_i)+1):length(x)]
    gametes <- gametes[nchar(gametes) > 0]
    gametes <- sites_matrix(gametes)
    out$gametes <- list(gametes)
  } else {
    out$positions <- list(numeric())
    out$gametes <- list(NULL)
  }

  if (has_tree(x))
    out$tree <- tree
  out
}

#' Parse MS output from R
#'
#' @param x character vector output from MS
#' @param include_seeds logical whether to include seeds
#'
#' Parse MS output from a character vector of lines of output to a tidy tibble.
#' Each replicate is a row, with a list-column for positions and a matrix of
#' gamete states.
#'
#' A \code{tree} column will be added if ms is run with \code{-T} (or, through 
#' command line arguments, \code{T=""}.
#'
#' You can use the purrr package to write map functions to summarize the
#' list-column of gametes.
#'
#' importFrom(magrittr,"%>%")
#' @export
#'
#' @examples
#' call_ms(10, 10, "-t 5") %>% parse_ms() 
parse_ms <- function(x, include_seeds=FALSE) {
  cmd <- x[1]
  seeds <- as.integer(strsplit(x[2], " ")[[1]])
  x <- x[4:length(x)]  # drop first few lines of metadata
  res_grp <- cumsum(x == "//")
  res <- split(x, res_grp)
  sims_lst  <- lapply(res, tidy_sim)
  sims <- do.call(rbind, sims_lst)
  if (include_seeds)
    sims$seeds <- paste(seeds, collapse=' ')
  # sims <- bind_rows(sims_lst)  # bind_rows() fails over 1000 entries
  sims$rep <- seq_along(sims_lst)

  colorder <- c('rep', 'segsites', 'positions', 'gametes')
  if (include_seeds)
    colorder <- c('rep', 'seeds', 'segsites', 'positions', 'gametes')
  if ('time' %in% colnames(sims)) {
    # unpack this special column
    sims <- mutate(sims, tmrca=map_dbl(time, first),
                         ttot=map_dbl(time, nth, n=2)) 
    colorder <- c(colorder, 'tmrca', 'ttot')
  }
  if ('tree' %in% colnames(sims))
    colorder <- c(colorder, 'tree')
  out <- sims[, colorder]
  out
}

