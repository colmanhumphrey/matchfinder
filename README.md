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

rows <- 600L
num_weight_vecs <- 25L
x_mat <- cbind(rnorm(rows),
               runif(rows))
treat_vec <- (1L:rows) %in%
    sample(1L:rows, size = floor(rows * 0.45),
           prob = (x_mat[, 2] + 2) / 5)
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

We're only generating 25 here, but for a serious go
of things, we'd probably want closer to 1,000.
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
    n_sinks = c(0L, 30L, 100L, 200L)
)

## 2. get all brier scores for all results
wr_briers <- lapply(all_wr_matches, function(by_sinks) {
    unlist(lapply(by_sinks, function(indiv_match_list) {
        brier_score_cv(
            x_mat = x_mat,
            match_list = indiv_match_list,
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
    n_sinks = c(0L, 30L, 100L, 200L),
    silent = TRUE
)
```

And indeed they're identical in terms of the matches:

```r
> ## to verify (can use any of the matches):
> all.equal(all_wr_matches[["100"]][[2]],
          brier_wr_matches[["matches_by_sinks"]][["100"]][[2]])
> [1] TRUE
```

The Brier scores however won't be the same, since they're based off
random draws to split training and testing. Note that by default, the
prediction method used is
[`xgboost`](https://cran.r-project.org/web/packages/xgboost/).  You
may want to change this. Both `brier_score` (and therefore
`brier_score_cv`) and `brier_[non]bipartite_matches` take
a parameter `match_predict_function`, which is by default
the result of calling
`match_predict_xgb()`: `match_predict_xgb` with all its defaults
(defaults for defaults!). This generates a function that only takes one
parameter, `train_test_list`; this is what `match_predict_function` is expected
to be. See `match_predict_xgb` and `match_predict_linear` for more details.
Both take parameters, so can be adjusted to your preference:

- if you're happy with `xgboost` but just want to change the
  parameters, you just call `match_predict_xgb` with your desired
  params:
  ```r
  some_result <- any_function_that_uses_match_predict_function(
      ...,
      ...,
      match_predict_function = match_predict_xgb(
          nrounds = 100,
          nthread = 4,
          params = list(eta = 0.3, max.depth = 6)
          ## anything else gets passed to xgboost::xgb.train
     ),
     ...
  )
  ```
  If you need more flexibility than that (e.g. you want to change
  `objective = "binary:logistic"` or anything), the functions produced by
  `match_predict_xgb`
  aren't too complex, you can build it manually
- if you want a simpler linear model - either logistic (`glm`) or actually
  linear (`lm`), pass `match_predict_linear(use_linear_lm = FALSE) # or true!`.
  It'll default to `glm`; to use `lm`, pass in `use_linear_lm = TRUE`.
- For a totally custom prediction function, again the internals of either should
  make clear what's happening

Anyway, now that we have our matches we want to choose one of them! For a given
number of sinks, we can default to using the hardest to predict match: the match
with the largest Brier score. But there are two things left to consider:
- For any given number of sinks, even if we choose the best, how good is this
  resulting match? Close to random, or still highly separable / easy to predict?
- If we have the best match per sink number, which match should we actually go
  with? What number of sinks is appropriate?

The corresponding paper goes into more detail, but essentially what we do is
score each match based on a permutation test. We could in theory just check if
the Brier score is above 0.25 or not, but this ignores many factors including
the overfitting behaviour of the chosen prediction function.

The permutation test:
- We take a given match
- For each iteration: we randomly choose the labelling for each pair (i.e.
  randomly choose which is treated vs control)
- We run the prediction algorithm and record the resulting Brier score
- We do this "some large number" of times
- The permutation score is how many permutation Briers were below the match's

This answers both questions:
- We can now effectively judge a match on a non-relative level
- We choose the number of sinks that first gives us an acceptable permutation
  score (going beyond that is throwing data away for higher variance in the
  actual result)

Let's now calculate this permutation score. In the below,
you could of course have used the two named elements of `brier_wr_matches`
either:

```r
perm_results <- permutation_matches(
    matches_by_sinks = all_wr_matches,
    briers_by_sinks = wr_briers,
    x_mat = x_mat,
    n_sinks = c(0L, 30L, 100L, 200L)
)
```

It's useful (but not necessary) to use `purrr` for the results, so let's load that up:

```r
library(purrr)
```

The resulting object has two named elements:

- `permutation_brier_scores`: this is a list with the same shape as
  `all_wr_matches` and `wr_briers`: that is, a list with an element per sink,
  and each of those is a result per match (in this case just a list of vectors). Let's look at this now (this is from a random run; you'll have different results):
  ```r
> perm_results$permutation_brier_scores
$`0`
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

$`30`
 [1] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 0.96 1.00 1.00 1.00 1.00 0.93 1.00 1.00
[16] 1.00 0.99 1.00 1.00 1.00 1.00 0.96 1.00 0.99 1.00

$`100`
 [1] 0.95 1.00 0.95 1.00 1.00 1.00 0.83 0.97 0.91 0.93 1.00 0.91 0.81 1.00 1.00
[16] 1.00 1.00 0.93 1.00 0.83 1.00 0.95 1.00 0.95 0.96

$`200`
 [1] 0.76 0.84 0.62 0.52 0.64 0.84 0.94 0.54 0.46 0.92 0.85 0.80 1.00 0.67 0.92
[16] 0.82 0.76 0.88 0.21 0.68 0.93 0.64 0.68 0.82 0.81
  ```
- `best_matches`: this gives a list with each element being a match, one per
  sink value. That is, we simply choose the best match by Brier score for each
  sink level and return it. Well that would be too little, we actually return a
  new list/object with the following elements:
  - `n_sinks` : somewhat simple here, try `map(perm_results$best_matches,
    "n_sinks")` or `map_int` if you want to get spooky
  - `raw_brier`: the raw Brier score from the input match. We already inputted
    the data that generated this result with `wr_briers`. We can see this with
    `lapply(wr_briers, max)` vs `map(perm_results$best_matches, "raw_brier")`
  - `permutation_brier`: the result of running the permutation Briers and counting
    how many scores were above the actual Brier. By default the process will only
    actually run the Brier scoring on the best outcome, **so we'll have
    the best raw Brier also gives the best permutation Brier, but this doesn't
    have to be true in general**. Hence why `lapply(wr_briers, which.max)`
    matches `lapply(perm_results$permutation_brier_scores, which.min)`.
  - `match_list`: the original match that we inputted into `permutation_matches`!
    Let's look at 30 for a second: `best_30_index <- which.max(wr_briers[["30"]])`.
    Then we'll see that `perm_results$best_matches[["30"]]$match_list` is identical
    to the input `all_wr_matches[["30"]][[best_30_index]]`.

Having the whole distribution of permutation scores available to us can be very
handy, in case we're worried we got a weird outlier etc. For example, in my case
I get the following results:

``` r
> map_dbl(perm_results$best_matches, "permutation_brier")
   0    30   100   200
1.00  0.93  0.81  0.21
```
In contrast, drawing a histogram can reveal what's really happening:

``` r
$`0`
0 |                             █| 1

$`30`
0 |                           ▁▁█| 1

$`100`
0 |                        ▂  ▂▃█| 1

$`200`
0 |      ▂      ▂ ▂▂ ▂▅█ ▅▂██▂█▂▂| 1
```

Here's the output for the four best matches:
``` r
    treat_est  treat_se num_matches num_sinks treat_sd raw_brier perm_brier
    0.50753   0.10678           270         0 1.7546   0.23732         1.00
    0.44120   0.10370           240        30 1.6066   0.24966         0.93
    0.34341   0.11819           170       100 1.5410   0.25219         0.81
    0.14635   0.15058            70       200 1.2598   0.27098         0.21
```

The simple regression, `y ~ treatment_vector`, gives an estimate of `0.57579`
with a standard error of `0.12375`. In this case, the full regression of `y ~
treatment_vector + X` is nearly optimal (but we only know that because we
generated the data), which gives an estimate of `0.43776` with a standard
error of `0.08504`.

So what do we see?

### Non Bipartite

```r
```
