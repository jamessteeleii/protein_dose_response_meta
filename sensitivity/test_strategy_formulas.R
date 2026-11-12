# Sampling-free unit checks for the protein-strategy formula transformation.

source("R/sensitivity_strategy.R")

find_target_command <- function(expression, target_name) {
  if (!is.call(expression)) return(NULL)
  if (
    identical(expression[[1L]], as.name("tar_target")) &&
      length(expression) >= 3L &&
      identical(expression[[2L]], as.name(target_name))
  ) {
    return(expression[[3L]])
  }
  for (element in as.list(expression)[-1L]) {
    result <- find_target_command(element, target_name)
    if (!is.null(result)) return(result)
  }
  NULL
}

pipeline <- parse(file = "_targets.R")
formula_command <- NULL
for (expression in pipeline) {
  formula_command <- find_target_command(
    expression, "model_formulas_primary"
  )
  if (!is.null(formula_command)) break
}
if (is.null(formula_command)) {
  stop("Could not locate the primary formula registry.", call. = FALSE)
}

primary_formulas <- eval(formula_command, envir = baseenv())
if (length(primary_formulas) != 26L) {
  stop("Expected 26 primary formulas.", call. = FALSE)
}

strategy_formulas <- vapply(seq_along(primary_formulas), function(index) {
  make_protein_strategy_formula(list(
    model_index = index,
    formula = primary_formulas[[index]]
  ))
}, character(1))

has_main <- grepl("\\bstrategy_code\\b", strategy_formulas)
has_dose_interaction <- grepl(
  "strategy_code:|by\\s*=\\s*strategy_code",
  strategy_formulas
)
expected_dose_interaction <- !seq_along(strategy_formulas) %in% c(1L, 5L)

stopifnot(
  all(has_main),
  identical(has_dose_interaction, expected_dose_interaction),
  all(vapply(strategy_formulas, function(x) {
    !inherits(try(stats::as.formula(x), silent = TRUE), "try-error")
  }, logical(1)))
)

# Model 26 already contains a protein-by-resistance-training smooth. The
# strategy sensitivity adds a separate strategy-difference smooth and must not
# inject a linear strategy term inside the existing RT smooth.
model_26_formula <- strategy_formulas[[26L]]
stopifnot(
  grepl(
    "s\\(protein_intake(?:_05)?,\\s*by\\s*=\\s*resistance_exercise_factor",
    model_26_formula,
    perl = TRUE
  ),
  grepl(
    "s\\(protein_intake(?:_05)?,\\s*by\\s*=\\s*strategy_code",
    model_26_formula,
    perl = TRUE
  ),
  !grepl(
    "s\\(protein_intake(?:_05)?\\s*\\+",
    model_26_formula,
    perl = TRUE
  ),
  !grepl("s\\([^)]*strategy_code:", model_26_formula)
)

# Exercise both formula-side remedies that are absent from the preregistered
# strings but can occur in the final 104-candidate registry.
remedy_cases <- list(
  list(
    model_index = 2L,
    formula = paste(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +",
      "protein_intake_05 + (1 | study_no) + (1 | arm_no) +",
      "(1 | outcome_no)"
    )
  ),
  list(
    model_index = 4L,
    formula = paste(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +",
      "splines::ns(protein_intake, df = 3) + (1 | study_no) +",
      "(1 | arm_no) + (1 | outcome_no)"
    )
  )
)
invisible(lapply(remedy_cases, function(specification) {
  stats::as.formula(make_protein_strategy_formula(specification))
}))

within_cases <- list(
  list(
    model_index = 1L,
    contrast_names = NULL,
    formula = paste(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +",
      "(1 | study_no) + (1 | arm_no) + (1 | outcome_no)"
    )
  ),
  list(
    model_index = 2L,
    contrast_names = "protein_contrast_1",
    formula = paste(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +",
      "protein_contrast_1 + (1 + protein_contrast_1 | study_no) +",
      "(1 | arm_no) + (1 | outcome_no)"
    )
  ),
  list(
    model_index = 3L,
    contrast_names = c("protein_contrast_1", "protein_contrast_2"),
    formula = paste(
      "yi|se(sqrt(vi)) ~ 0 + Intercept + resistance_exercise_code +",
      "protein_contrast_1 + protein_contrast_2 +",
      "protein_contrast_1:resistance_exercise_code +",
      "protein_contrast_2:resistance_exercise_code +",
      "(1 + protein_contrast_1 + protein_contrast_2 | study_no) +",
      "(1 | arm_no) + (1 | outcome_no)"
    )
  )
)
invisible(lapply(within_cases, function(case) {
  transformed <- make_protein_strategy_within_formula(
    formula = case$formula,
    model_index = case$model_index,
    contrast_names = case$contrast_names
  )
  stopifnot(grepl("strategy_code", transformed, fixed = TRUE))
  if (!case$model_index %in% c(1L, 5L)) {
    stopifnot(all(vapply(case$contrast_names, function(name) {
      grepl(paste0(name, ":strategy_code"), transformed, fixed = TRUE)
    }, logical(1))))
  }
  stats::as.formula(transformed)
}))

cat(
  paste(
    "PASS: 26 arm-based formulas, half-unit/natural-spline remedies,",
    "and within-contrast strategy formulas.\n"
  )
)
