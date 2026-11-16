source("R/prior_sensitivity.R")

# For a standard Normal prior and likelihood N(1, 1), the posterior
# under prior power a is N(1 / (1 + a), 1 / (1 + a)).
set.seed(101)
n <- 24000L
chain <- rep(1:4, each = n / 4L)
draw <- stats::rnorm(n, mean = 0.5, sd = sqrt(0.5))
log_prior <- stats::dnorm(draw, log = TRUE)
original <- prior_sensitivity_weights(log_prior, chain, 1)
stopifnot(original$reliable, all(original$weights == 1 / n),
          isTRUE(all.equal(sum(original$weights * draw), mean(draw))))
for (power in c(0.8, 1.2)) {
  ps <- prior_sensitivity_weights(log_prior, chain, power)
  stopifnot(ps$reliable, abs(sum(ps$weights) - 1) < 1e-12,
            abs(sum(draw * ps$weights) - 1 / (1 + power)) < 0.015)
}

# Strong reweighting with inadequate overlap should fail the reliability gate.
unstable <- suppressWarnings(prior_sensitivity_weights(
  c(log_prior[-1L], -500), chain, 0.8
))
stopifnot(!unstable$reliable)

# Both endpoints of the equivalence region are inclusive.
result <- summarise_prior_sensitivity_contrast(
  c(-0.11, -0.10, 0, 0.10, 0.11), rep(0.2, 5)
)
stopifnot(isTRUE(all.equal(result$equivalence, 0.6)),
          result$gain == 0.2, result$loss == 0.2,
          result$conclusion == "uncertain", result$low == -0.11, result$high == 0.11)
result <- summarise_prior_sensitivity_contrast(c(0, 0.11), c(0.95, 0.05))
stopifnot(result$conclusion == "equivalence")

cat("Prior-sensitivity tests passed.\n")
