# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

iterate_cyclic <- function(X, P, states, threshold) {
    .Call(`_calib_iterate_cyclic`, X, P, states, threshold)
}

iterate_index <- function(X, P, R) {
    .Call(`_calib_iterate_index`, X, P, R)
}

iterate_life <- function(X, P, rule) {
    .Call(`_calib_iterate_life`, X, P, rule)
}

iterate_total <- function(X, P, R) {
    .Call(`_calib_iterate_total`, X, P, R)
}

