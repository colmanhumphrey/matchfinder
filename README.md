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

### Bipartite

First we'll load the library and setup some data:
```r
library(matchfinder)

treat_effect <- 0.3

rows <- 500L
num_weight_vecs <- 5L
x_mat <- cbind(rnorm(rows),
               runif(rows))
treat_vec <- (1L:rows) %in% fixed_sample(1L:rows, floor(rows * 0.45))
y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * treat_effect + rnorm(rows)
```

Next let's compute the covariance matrix (and distance matrix).
Note that `covariance_with_ranks` allows us to specify
which columns should be used on a rank level if desired.
```r
cov_x <- covariance_with_ranks(x_mat)

## we don't use this, just showing how it's done
## under the hood within the higher-level functions
## or if you wanted to call `bipartite_matches` directly
dist_mat <- weighted_mahal(x_mat,
    cov_x = cov_x,
    weight_vec = c(0.66, 0.33),
    treat_vec = treat_vec
)
```

We generate some random weight vectors:
idea is that we'll end up using the one that
gave us the best match.
You don't have to use this function, it just
has some reasonble features for generating random
vectors, including setting minimum weights, prior
weights etc.
```r
weight_vecs <- generate_random_weights(
    prior_weights = c(2, 1),
    number_vectors = num_weight_vecs,
    minimum_weights = c(0.1, 0.1)
)
```

We'll now generate many matches. Here we only try two different
sink values; in general we can either try many more here,
or increase until we get the desired randomness in the match.

We'll do this two ways. Firstly, the two steps explicitly:

1. Generate the matches
2. Compute the Brier score for each match
```r
## 1. Here we generate all the matches: one for each pair
##    of weight vector and sink value
all_wr_matches <- all_bipartite_matches(
    x_mat = x_mat,
    cov_x = cov_x,
    weight_list = weight_vecs,
    treat_vec = treat_vec,
    match_method = "with_replacement",
    n_sinks = c(0L, 30L)
)

## 2. get all brier scores for all results
wr_briers <- lapply(all_wr_matches, function(by_sinks) {
    unlist(lapply(by_sinks, function(indiv_match_list) {
        brier_score_cv(
            x_mat = x_mat,
            match_list = indiv_match_list
        )
    }))
})
```

Or we can call `brier_bipartite_matches`, which performs both steps for us:
```r
brier_wr_matches <- brier_bipartite_matches(
    x_mat = x_mat,
    cov_x = cov_x,
    weight_list = weight_vecs,
    treat_vec = treat_vec,
    match_method = "with_replacement",
    n_sinks = c(0L, 30L),
    silent = TRUE
)

## to verify (can use any of the matches):
all.equal(all_wr_matches[["0"]][[2]],
          brier_wr_matches[["matches_by_sinks"]][["0"]][[2]])
```


### Non Bipartite

```r
```
