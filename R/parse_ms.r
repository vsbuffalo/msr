
#' Call MS from R
#'
#' Call Hudon's MS from R, returning a string of results (that can be parsed
#' with parse_ms()). 
#'
#' @param nsam number of samples (gametes) to draw
#' @param howmany how many replicates to run
#' @param cmd the command to pass to MS 
#' 
#' @export
call_ms <- function(nsam, howmany, cmd, verbose=TRUE) {
  ms_cmd <- sprintf("ms %d %d %s", nsam, howmany, cmd)
  if (verbose)
    message(sprintf("command: %s", ms_cmd))
  system(ms_cmd, intern=TRUE)
}


#' A "tee" (in Unix sense) that writes MS's results to file
#'
#' @param x MS results @param con to save output to
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
write_tee <- function(x, con) {
  writeLines(x, con=con)
  return(x)
}

#' Parse MS's key/value pairs, e.g. segsites and positions
#' returning the values as list if there are more than one
#' @keywords internal
parse_keyvals  <- function(x) {
  keyvals <- gsub("(\\w+): +(.*)", "\\1;;\\2", x, perl=TRUE)
  tmp <- strsplit(keyvals, ";;")[[1]]
  key <- tmp[1]
  vals <- as.numeric(strsplit(tmp[2], " ")[[1]])
  if (length(vals) == 1)
    return(vals)
  else
    return(list(vals))
}


#' Convert gaemtes' alleles intro matrix
#' @keywords internal
sites_matrix <- function(x) {
  do.call(rbind, lapply(x, function(y) as.integer(unlist(strsplit(y, "")))))
}

#' Tidy a single simulation result from MS
#' @keywords internal
tidy_sim <- function(x) {
  # first element is the delimiter "//", last element is blank line (except for
  # last sim) 
  segsites <- parse_keyvals (x[2])
  positions <- parse_keyvals (x[3])
  gametes <- x[4:length(x)]

  # remove empty line if there
  gametes <- gametes[nchar(gametes) > 0]
  gametes <- sites_matrix(gametes)
  tibble(segsites, positions=positions, gametes=list(gametes))
}

#' Parse MS output from R
#'
#' @param x character vector output from MS
#'
#' Parse MS output from a character vector of lines of output to a tidy tibble.
#' Each replicate is a row, with a list-column for positions and a matrix of
#' gamete states.
#'
#'
#' You can use the purrr package to write map functions to summarize the
#' list-column of gametes.
#'
#' importFrom(magrittr,"%>%")
#' @export
#'
#' @examples
#' call_ms(10, 10, "-t 5") %>% parse_ms() 
parse_ms <- function(x) {
  cmd <- x[1]
  seeds <- as.integer(strsplit(x[2], " ")[[1]])
  x <- x[4:length(x)]  # drop first few lines of metadata
  res_grp <- cumsum(x == "//")
  res <- split(x, res_grp)
  sims_lst  <- lapply(res, tidy_sim)
  sims <- do.call(rbind, sims_lst)
  # sims <- bind_rows(sims_lst)  # bind_rows() fails over 1000 entries
  sims$rep <- seq_along(sims_lst)
  sims[, c('rep', 'segsites', 'positions', 'gametes')]
}

