# Structure

1. We have `weighted_mahal`, to get the weighted mahalanobis distances (which relies on `covariance_with_ranks`). Note that this can be rank-adjusted. In file `mahal.R`.



N. all_weighted_mahal takes a weight_list of weighted vecs:
FOR EACH WEIGHT VEC in weight_list
 - `mahal - weighted_mahal` to get the distances
 - OPTIONALLY adds a caliper using `caliper - add_caliper`
 - cases:
   - **bipartite, with replacement**: 
   - **bipartite, no replacement**:
   - **nonbipartite**, with/without replacement: both call `nonbipartite_match`
