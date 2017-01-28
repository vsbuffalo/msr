
#' Calculate Pi, via the sum of site heterozygosities 
#' 
#' @param x a gamete matrix
#' @export
theta_pi <- function(x) { p <- colSums(x)/nrow(x); sum(2*p*(1-p)) } 

#' Calculate Watterson's Theta
#' 
#' @param s number of segregating sites
#' @param n number of samples
#' @export
theta_W <- function(s, n) s/(sum(1/(1:(n-1))))

tajD_num <- function(pi, s, n) {
  a1 <- sum(1/(1:(n-1)))
  pi-s/a1
} 

tajD_denom <- function(s, n) {
  a1 <- sum(1/(1:(n-1)))
  a2 <- sum(1/(1:(n-1))^2)
  b1 <- (n+1)/(3*(n-1))
  b2 <- (2*(n^2 + n + 3))/(9*n*(n-1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n+2)/(a1*n) + a2 / a1^2
  e1 <- c1/a1 
  e2 <- c2/(a1^2 + a2)
  (sqrt(e1*s) + e2*s*(s-1))
}


#' Calculate Tajima's D
#' 
#' @param pi theta pi, or pairwise differences
#' @param s numbr of segregating sites
#' @param n number of samples
#' @export
tajD <- function(pi, s, n) {
  a1 <- sum(1/(1:(n-1)))
  a2 <- sum(1/(1:(n-1))^2)
  b1 <- (n+1)/(3*(n-1))
  b2 <- (2*(n^2 + n + 3))/(9*n*(n-1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n+2)/(a1*n) + a2 / a1^2
  e1 <- c1/a1 
  e2 <- c2/(a1^2 + a2)
  D <- (pi-s/a1)/(sqrt(e1*s) + e2*s*(s-1))
  D
}


#' Mutate a tibble of results from parse_ms(), adding summary statistics
#' 
#' @param x tibble of results from parse_ms()
#' @param .n number of samples
#'
#' @export
sample_stats <- function(x, .n) {
  x %>% mutate(theta_pi=map_dbl(gametes, theta_pi), theta_W=theta_W(segsites, .n), D=tajD(theta_pi, segsites, .n))
               # D_num=tajD_num(theta_pi, segsites, .n), D_denom=tajD_denom(segsites, .n))
}


