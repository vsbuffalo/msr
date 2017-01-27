# Call/Parse Hudson's MS coalescent simulator results from R

This is a tiny package to call/parse MS results from R. I got sort of sick of
saving MS output to file and loading it into R each time I wanted a quick
graphic. The design is in following with Hadley's
[purrr](https://github.com/hadley/purrr) and
[dplyr](https://github.com/hadley/dplyr)-driven tidy way of manipulating
dataframes (or in this case, tibbles). 

## Example: Calculating Pi

Here's a simple example:

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
> res %>% summarize(theta_pi=mean(pi), theta_W=mean(segsites)/sum(1/(1:9)))
# A tibble: 1 × 1
      pi
   <dbl>
1 5.0866
```

The list-column approach stores the site matrices for each run in a `gametes`
list-column. Each element of the list is a matrix. Similarly, `positions` is a
list-column, where each element is a vector of positions. Calculations on the
gamete site matrix can be done using the `mutate(stat = map_dbl(gametes,
some_fun))` dplyr/purrr pattern (see the example above).

## Writing Output to File

It's likely for reproducibility's sake you'd like to save the MS run to file.
You can do this, still sending the output to the pipe with `write_tee()`:

```{R}
> call_ms(10, 100, "-t 5") %>% write_tee("ms.out") %>%
   parse_ms() %>% mutate(pi=map_dbl(gametes, pi))

# A tibble: 100 × 5
     rep segsites  positions         gametes    pi
   <int>    <dbl>     <list>          <list> <dbl>
1      1        8  <dbl [8]>  <int [10 × 8]>  1.86
2      2       13 <dbl [13]> <int [10 × 13]>  2.58
3      3        4  <dbl [4]>  <int [10 × 4]>  1.64
4      4       14 <dbl [14]> <int [10 × 14]>  3.86
5      5       20 <dbl [20]> <int [10 × 20]>  4.98
6      6        9  <dbl [9]>  <int [10 × 9]>  3.10
7      7        3  <dbl [3]>  <int [10 × 3]>  1.08
8      8       11 <dbl [11]> <int [10 × 11]>  3.10
9      9       23 <dbl [23]> <int [10 × 23]> 10.32
10    10       28 <dbl [28]> <int [10 × 28]> 10.02
# ... with 90 more rows
```
