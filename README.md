# matchfinder

<!-- badges: start -->
<!-- badges: end -->

`matchfinder` is for finding best matches based on predictability. The
general idea: a given match is bad if we can predict which units are
treated and which are control.

## Main Concept

First let's discuss the case of matched pairs:

- We have a set of `T` treatment units and `C` control units,
  for a total of `N = T + C` units
- We form `n` pairs, where `n ≤ T, n ≤ C`
  - Each pair consists of one treated unit and one control unit
  - Typically we restrict each treated unit to show up at most once in all pairs
  - Optionally we can restrict the control units the same way
- We may leave out some or many of the original units
- The groups of units, all resulting matched treated units, and all resulting
  matched control units, shouldn't differ systematically, to reduce bias
- The units within each pair should be similar, to reduce variance

However, it isn't enough to evaluate group-level statistics for the
bias part: say we have some variable `x` and in our pairs of (treated,
control), it only takes the values (0, 1), (1, 2), (2, 3) and (3, 0),
each 25% of the time. Then the variable will have the same group-level
statistics for treated and control units, but there's a clear pattern
in the data. If `x` has some non-linear effect on the outcome, this
would bias our results.

If you were given half the pairs, along with all relevant info about
the pairs (values of `x`, which unit is treated etc), you could
quickly learn the pattern given and easily predict which unit is the
treated when looking at the other half. This means it's easy to predict
which units are the treated units, and our match is bad.

## Installation

We're not yet on [CRAN](https://CRAN.R-project.org), but you can install
from github with `devtools` or `remotes`:

``` r
## install.packages("devtools")
devtools::install_github("ColmanHumphrey/matchfinder")

## OR:

## install.packages("remotes")
remotes::install_github("ColmanHumphrey/matchfinder")
```

## Example

We'll work through two small examples, one bipartite, one
non-bipartite.

```r
## We'll load the library in either case
library(matchfinder)
```

### Bipartite

```r
    rows <- 100L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             45L)
    y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * 0.3
    cov_x <- covariance_with_ranks(x_mat)

    some_dist_mat <- weighted_mahal(x_mat,
                                    cov_x = cov_x,
                                    weight_vec = c(0.66, 0.33),
                                    treat_vec = treat_vec)

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    all_wr_matches <- all_bipartite_matches(x_mat = x_mat,
                                            cov_x = cov_x,
                                            weight_list = weight_vecs,
                                            treat_vec = treat_vec,
                                            match_method = "with_replacement",
                                            n_sinks = c(0L, 4L))

    sink_brier_wr_matches <- sink_brier_bipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        treat_vec = treat_vec,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L),
        silent = TRUE)

```


### Non Bipartite
