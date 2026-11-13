# Synthetic tests only: no targets execution, fitted objects, compilation or sampling.
stopifnot(getRversion() == package_version("4.4.3"))
suppressPackageStartupMessages(library(tidyverse))
for (f in c("R/functions.R", "R/primary_model_averaging.R", "R/verification_artifacts.R",
  "R/age_duration_interactions.R", "R/age_duration_recovery_reporting.R",
  "R/age_duration_batch.R", "R/age_duration_family_targets.R")) source(f)
expect_error <- function(x) stopifnot(inherits(tryCatch({force(x); NULL}, error = identity), "error"))

# Historical spline frames keep the study factor instead of study_no. Use its
# labels while still rejecting other studies, row permutations or likelihoods.
identity_reference <- tibble::tibble(study_no = c("10", "20", "10"), arm_no = 1:3,
  outcome_no = 1:3, yi = c(.1, .2, .3), vi = c(.01, .02, .03))
legacy_frame <- dplyr::select(identity_reference, -study_no)
legacy_frame$study_factor <- factor(identity_reference$study_no, levels = c("20", "10"))
legacy_before <- legacy_frame
stopifnot(iad_same_observations(legacy_frame, identity_reference),
  iad_same_observations(identity_reference, legacy_frame), identical(legacy_frame, legacy_before),
  !iad_same_observations(legacy_frame[3:1, ], identity_reference))
for (column in c("study_factor", "arm_no", "outcome_no", "yi", "vi")) {
  changed <- legacy_frame
  if (column == "study_factor") changed$study_factor[1L] <- "20" else
    changed[[column]][1L] <- changed[[column]][1L] + 1
  stopifnot(!iad_same_observations(changed, identity_reference))
}
stopifnot(!iad_same_observations(dplyr::select(legacy_frame, -study_factor), identity_reference))

data <- tibble::tibble(study_no = "1", arm_no = 1, outcome_no = 1, time1 = 0, time2 = 12,
  yi = 0.1, vi = 0.1, protein_intake = 0, protein_intake_raw = 100,
  age_centre = 0, duration_centre = 0, resistance_exercise_code = 0)
fixture <- list(); specifications <- list()
for (family in iad_families()) {
  parents <- if (family == "relative") c(6L,7L,8L,12L,13L,14L,18L,19L,20L,24L,25L,26L) else
    c(6L,8L,12L,14L,20L,24L,25L,26L)
  fixture[[family]] <- list(data = stats::setNames(rep(list(data),4L),iad_outcomes()),
    parents = tidyr::crossing(outcome = iad_outcomes(),model_index = parents) |>
      dplyr::mutate(expected_bytes = 1024^3))
  for (outcome in iad_outcomes()) {
    index <- 26L
    for (interaction in c("age","duration","both")) for (parent in parents) {
      index <- index + 1L
      specifications[[length(specifications)+1L]] <- list(family = family,outcome = outcome,
        model_index = index,parent_index = parent,interaction = interaction,stage = "initial")
    }
  }
}
stopifnot(length(specifications)==336L)
pilots <- dplyr::bind_rows(iad_pilot_specs(specifications)) |>
  dplyr::mutate(diagnostic_pass=TRUE,data_sha256=iad_fingerprint(data))
for(family in iad_families()) iad_check_family_pilot(pilots,family,fixture,specifications)
failed <- pilots; failed$diagnostic_pass[failed$family=="relative"][1L] <- FALSE
expect_error(iad_check_family_pilot(failed,"relative",fixture,specifications))
iad_check_family_pilot(failed,"total",fixture,specifications)
iad_check_family_pilot(failed,"supplement",fixture,specifications)
expect_error(iad_check_family_pilot(pilots[-1L,],"relative",fixture,specifications))
expect_error(iad_check_family_pilot(dplyr::bind_rows(pilots,pilots[1L,]),"relative",fixture,specifications))
wrong <- pilots; wrong$data_sha256[wrong$family=="total"] <- "wrong"
expect_error(iad_check_family_pilot(wrong,"total",fixture,specifications))

Sys.setenv(PROTEIN_INTERACTION_PHASE="batch_fit")
stopifnot(iad_phase("pilot"),iad_phase("fit"),!iad_phase("loo"),!iad_phase("remedy"),!iad_phase("report"))
iad_cue <- function(phase) targets::tar_cue(mode=if(iad_phase(phase))"thorough"else"never")
graph <- unlist(lapply(iad_families(),iad_family_targets),recursive=FALSE)
stopifnot(!anyDuplicated(vapply(graph,function(x)x$name,character(1))),
  all(vapply(graph,function(x)identical(x$settings$error,"continue"),logical(1))))
for(family in iad_families()) {
  nodes <- stats::setNames(graph,vapply(graph,function(x)x$name,character(1)))
  fit <- nodes[[paste0("iad_",family,"_fit_file")]]
  deps <- all.names(fit$command$expr)
  stopifnot(paste0("iad_",family,"_pilot_gate_file") %in% deps,
    !any(paste0("iad_",setdiff(iad_families(),family),"_pilot_gate_file") %in% deps))
}
Sys.setenv(PROTEIN_INTERACTION_PHASE="none")

# Credit only actual saved initial fits when resuming; no fitting is performed.
test_storage <- function() {
  root <- tempfile("iad-batch-storage-")
  original <- iad_path
  assign("iad_path",function(...,archive=FALSE)file.path(root,...),envir=.GlobalEnv)
  on.exit(assign("iad_path",original,envir=.GlobalEnv),add=TRUE)
  plan <- iad_batch_storage_plan(specifications,fixture,iad_families(),600)
  stopifnot(all(plan$storage_pass),sum(plan$candidates)==336L,
    unique(plan$total_required_gib)==524)
  low <- iad_batch_storage_plan(specifications,fixture,iad_families(),500)
  stopifnot(!any(low$storage_pass))
  first <- specifications[[1L]]
  path <- iad_path(first$family,first$outcome,first$stage,sprintf("model_%02d.qs",first$model_index),archive=TRUE)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE); writeLines("fixture",path)
  resumed <- iad_batch_storage_plan(specifications,fixture,iad_families(),600)
  stopifnot(sum(resumed$saved)==1L,unique(resumed$total_required_gib)==522.5)
  expect_error(iad_batch_storage_plan(specifications,fixture,"unknown",600))
}
test_storage()

# Regression for the recovery reporting failure at modifyList(settings, ...).
grid <- list(keys=tibble::tibble(resistance_exercise_code=c(0,0,1,1),
  protein_intake=c(50,100,50,100)),grid_step=50,protein_min=50)
record <- list(outcome="whole_body",models=tibble::tibble(model_name="fixture",stacking_weight=1))
draws <- matrix(seq_len(40)/100,nrow=10)
average <- iad_recovery_model_average(record,grid,draws)
stopifnot(is.list(average$settings),average$settings$n_draws==10L)
support <- tibble::tibble(resistance_exercise_code=c(0,0,1,1),protein_intake_raw=c(50,100,50,100))
contrasts <- compute_primary_heatmap_contrasts(average,support,protein_min=50)
stopifnot(nrow(contrasts$contrasts)==8L)

# Test the destructive-step ordering with stubs; no real pruning or targets.
test_launcher_order <- function() {
  environment <- new.env(parent=.GlobalEnv)
  sys.source("sensitivity/run_age_duration_batch.R",envir=environment)
  events <- character()
  reject_second_plan <- FALSE
  environment$source <- function(file,local=FALSE,...) {
    caller <- parent.frame()
    assign("iad_storage_plan",function(key,verify=TRUE) {
      events <<- c(events,paste0("verify:",key))
      if(reject_second_plan && key=="total_strategy") stop("fixture invalid second plan")
      list(summary=data.frame(key=key,reclaimable_gib=1))
    },envir=caller)
    assign("iad_prune_completed_family",function(key,execute=FALSE) {
      stopifnot(isTRUE(execute));events <<- c(events,paste0("prune:",key))
    },envir=caller)
  }
  environment$run_age_duration_checkpoint <- function(stage,families,workers) {
    stopifnot(stage=="initial_fit",workers==4L,setequal(families,iad_families()))
    events <<- c(events,"initial_fit")
  }
  environment$run_age_duration_initial_batch(4L,TRUE)
  stopifnot(identical(events,c("verify:supplement_strategy","verify:total_strategy",
    "prune:supplement_strategy","prune:total_strategy","initial_fit")))
  events <- character();reject_second_plan <- TRUE
  expect_error(environment$run_age_duration_initial_batch(4L,TRUE))
  stopifnot(!any(grepl("prune:|initial_fit",events)))
  events <- character()
  environment$run_age_duration_initial_batch(4L,FALSE)
  stopifnot(identical(events,"initial_fit"))
}
test_launcher_order()

test_clean_preparation <- function() {
  runner <- new.env(parent = .GlobalEnv)
  sys.source("sensitivity/run_age_duration_batch.R", envir = runner)
  has_recovery <- FALSE
  runner$file.exists <- function(path) {
    stopifnot(basename(path) == "accepted_manifest.csv")
    has_recovery
  }
  script <- body(runner$run_age_duration_checkpoint)
  # Intercept the local stage runner before it can call targets.
  runner$stages_run <- list()
  runner$record_stage <- function(targets, phase) {
    runner$stages_run <- c(runner$stages_run, list(c(phase, targets)))
  }
  define_runner <- which(vapply(as.list(script), function(expr) {
    is.call(expr) && identical(expr[[1L]], quote(`<-`)) &&
      identical(expr[[2L]], quote(run))
  }, logical(1L)))
  stopifnot(length(define_runner) == 1L)
  script[[define_runner]] <- quote(run <- record_stage)
  body(runner$run_age_duration_checkpoint) <- script
  preparation <- c("none", "iad_preparation_files", "iad_batch_storage_gate_file")
  runner$run_age_duration_checkpoint("prepare")
  stopifnot(identical(runner$stages_run, list(preparation)))
  runner$stages_run <- list()
  has_recovery <- TRUE
  runner$run_age_duration_checkpoint("prepare")
  stopifnot(identical(runner$stages_run, list(
    c("recover_report", "iad_recovery_reporting_audit_files"), preparation
  )))
}
test_clean_preparation()
cat("Independent pilot gates, batch phase boundaries, storage/resume guards and recovery reporting regression passed. No pipeline run.\n")
