# Call/Parse Hudson's MS coalescent simulator results from R

This is a tiny package to call/parse MS results from R. I got sort of sick of
saving output to file when the end result was a quick graphic. The design is in
following with Hadley's [purrr](https://github.com/hadley/purrr) and
[dplyr](https://github.com/hadley/dplyr)-driven tidy way of manipulating
dataframes (or in this case, tibbles). Here's a simple example:

```{R}
> library(msr)
> library(tidyverse)
> pi <- function(x) { p <- colSums(x)/nrow(x); sum(2*p*(1-p)) } 
> res <- call_ms(10, 100, "-t 5") %>% parse_ms() %>% mutate(pi=map_dbl(gametes, pi))
> res
# A tibble: 100 × 5
     rep segsites  positions         gametes    pi
   <int>    <dbl>     <list>          <list> <dbl>
1      1       19 <dbl [19]> <int [10 × 19]>  6.06
2      2       14 <dbl [14]> <int [10 × 14]>  3.64
3      3       21 <dbl [21]> <int [10 × 21]>  6.58
4      4       21 <dbl [21]> <int [10 × 21]>  6.20
5      5       12 <dbl [12]> <int [10 × 12]>  4.56
6      6        4  <dbl [4]>  <int [10 × 4]>  1.56
7      7       13 <dbl [13]> <int [10 × 13]>  4.10
8      8        9  <dbl [9]>  <int [10 × 9]>  2.62
9      9       16 <dbl [16]> <int [10 × 16]>  5.92
10    10       12 <dbl [12]> <int [10 × 12]>  3.32
# ... with 90 more rows
> res %>% summarize(pi=mean(pi))
# A tibble: 1 × 1
      pi
   <dbl>
1 5.0866
```

