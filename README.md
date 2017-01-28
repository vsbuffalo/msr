# Call/Parse/Analyze Hudson's MS coalescent simulator results from R

This is a tiny package to call/parse MS results from R. I got sort of sick of
saving MS output to file and loading it into R each time I wanted a quick
graphic. The design is in following with Hadley's
[purrr](https://github.com/hadley/purrr) and
[dplyr](https://github.com/hadley/dplyr)-driven tidy way of manipulating
dataframes (or in this case, tibbles). 

## Example: Calculating Theta Pi and Watterson's Theta

Here's a simple example:

```{R}
> library(msr)
> library(tidyverse)
> res <- call_ms(10, 100, t=5) %>% parse_ms() %>% mutate(pi=map_dbl(gametes, theta_pi))
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
```

Note that `call_ms()` (and its sibling function `ms()`, which also
automatically parses the results) take arguments in two ways: (1) as a command
line string, and (2) as function arguments (as theta was specified above). For
the latter, multiple-valued arguments are passed as a vector, e.g. to run ms
with a theta of 30 (`-t 30`), a rho of 10 for 1000 sites (`-r 10 1000`),
sampling 10 gametes, and replicated 100 times use: `ms(10, 1000, t=30, r=c(10,
1000))`.  Taking arguments directly through the function like this makes
calling MS with multiple parameters easier, through the purrr function
`invoke_rows()` (see further below).

Returning to our example above, we can now analyze these results, using the
packaged summary statistics functions `theta_pi()` and `theta_W()` on the data
(though see below for an alternate way with the function `sample_stats()`:

```{R}
> res <-  res %>% mutate(theta_pi=map_dbl(gametes, theta_pi), theta_W=theta_W(segsites, 10))
> res 
# A tibble: 100 × 7
     rep segsites  positions         gametes        pi  theta_pi   theta_W
   <int>    <dbl>     <list>          <list>     <dbl>     <dbl>     <dbl>
1      1       21 <dbl [21]> <int [15 × 21]>  6.506667  6.506667  7.423201
2      2       14 <dbl [14]> <int [15 × 14]>  4.817778  4.817778  4.948801
3      3       11 <dbl [11]> <int [15 × 11]>  3.537778  3.537778  3.888343
4      4       12 <dbl [12]> <int [15 × 12]>  2.293333  2.293333  4.241829
5      5       24 <dbl [24]> <int [15 × 24]>  7.431111  7.431111  8.483658
6      6       44 <dbl [44]> <int [15 × 44]> 16.551111 16.551111 15.553374
7      7       20 <dbl [20]> <int [15 × 20]>  5.066667  5.066667  7.069715
8      8       15 <dbl [15]> <int [15 × 15]>  2.275556  2.275556  5.302286
9      9       18 <dbl [18]> <int [15 × 18]>  5.617778  5.617778  6.362744
10    10       19 <dbl [19]> <int [15 × 19]>  7.253333  7.253333  6.716229
# ... with 90 more rows
```

Now, we find the mean of these summary statistics:

```{R}
> res %>% summarize(theta_pi = mean(theta_pi), theta_W = mean(theta_W))
 A tibble: 1 × 2
  theta_pi  theta_W
     <dbl>    <dbl>
1 4.950222 6.062281
```

The list-column approach stores the site matrices for each run in a `gametes`
list-column. Each element of the list is a matrix. Similarly, `positions` is a
list-column, where each element is a vector of positions. Calculations on the
gamete site matrix can be done using the `mutate(stat = map_dbl(gametes,
some_fun))` dplyr/purrr pattern (see the example above).

## Sample Statistics

Some basic sample statistics functions like those in MS's `sample_stats`
package are included in the package: `theta_W()`, `theta_pi()`, `tajD()`, and a
function `sample_stats()` that calculates these in a pipeline. 

To make quick runs simple, `ms()` runs both `call_ms()` and `parse_ms()`. We
combine this with `sample_stats()` below:

```{R}
ms(nsam=10, howmany=50, t=30) %>% sample_stats(.n=10) %>% 
    summarize(theta_pi = mean(theta_pi), theta_W = mean(theta_W))
```

Or, a quick graphic example:

```{R}
ggplot(ms(10, 1000, t=30) %>% sample_stats(.n=10)) + geom_histogram(aes(x=D))
```


## Writing Output to File

It's likely for reproducibility's sake you'd like to save the MS run to file.
You can do this, still sending the output to the pipe with `write_tee()`:

```{R}
> call_ms(10, 100, "-t 5") %>% write_tee("ms.out") %>%
   parse_ms() %>% mutate(pi=map_dbl(gametes, theta_pi))

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

## MS Runs Across Multiple Parameters

The purrr package has a really nice feature: if we have rows of a dataframe
representing a set of parameters, we can "invoke" a function on these
parameters with the function `invoke_row()` (or `by_rows()`). With `msr`, this
is useful if you want to analyze the results of multiple MS runs across
different parameter values. While this is possible with other approaches in R,
it's much easier using the tidy data approach. Below, I replicate Table 1 of
Kevin Thornton's neat paper [*Recombination and the Properties of Tajima's D in
the Context of Approximate-Likelihood
Calculation*](http://www.genetics.org/content/171/4/2143) using this approach
(though do note that I am using N=10^4 rather than N=10^5 as in the paper since
our department's happy hour is soon, and I didn't want to wait for simulations
to complete for this example):

```{R}
> library(tidyverse)
> library(msr)

# generate the parameter tibble:
> params <- tibble(nsam=30, howmany=10^4, t=20, r=list(c(0, 1e3), c(10, 1e3), c(50, 1e3))) 

# run the simulations
> res <- params %>% invoke_rows(.f=ms, .to="msout")

# append a sample_stats list-column, calculating sample_stats on each MS run,
# then unnest, bringing these columns into the main tibble
> stats <- res %>% mutate(rho=map_dbl(r, first)) %>% 
              mutate(sample_stats=map(msout, sample_stats, .n=first(nsam))) %>% unnest(sample_stats)

# group by all changing parameters, and calc summaries of the summary statistics
> stats %>% group_by(rho) %>% summarize(ED=mean(D), sd_D=sd(D))
# A tibble: 3 × 3
    rho          ED       sd_D
  <dbl>       <dbl>      <dbl>
1     0 -0.06608313 0.18045349
2    10 -0.04190666 0.13234880
3    50 -0.02816332 0.08536906
```

The second column `ED` is reasonably close to Thornton's second column of Table
1 in his paper.
