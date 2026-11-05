##### scr to run components for hybrid vr study #####
# created by cas 28 jul 2026 #
# modified by cas 22 oct 2026
#  COPdiff-LFdiff correlation analysis: real-world (IRL) and hybrid-VR (VR)
#
# Excel layout expected for each workbook:
#   - Row 1 = participant ID
#   - Rows 2-26 =  participant's 25 pre-switch trials
#   - The five blocks in order with five measurements per block
#
# Required files in data_dir:
#   pre_5_cop_vr.xlsx, pre_5_cop_irl.xlsx
#   pre_5_lf_vr.xlsx,  pre_5_lf_irl.xlsx

data_dir <- "/Users/motorskill/Desktop/R01s/vr/rep_vr/data_analyses"

#####
# Choose the analysis window by setting window_name to ONE option 
window_name <- "last3_each_block"

# Principal analysis: final three lifts in each of five blocks (15 values):
trial_indices <- c(3:5, 8:10, 13:15, 18:20, 23:25)

# secondary analysis: all five lifts from the final three blocks (15 values):
# window_name <- "last3_blocks"
# trial_indices <- 11:25

# secondary analysis: all 25 pre-switch measurements:
# window_name <- "all25"
# trial_indices <- 1:25

# To switch windows, comment out the active window_name/trial_indices lines above
#####

files <- list(
  VR = list(cop = file.path(data_dir, "pre_5_cop_vr.xlsx"),
            lf  = file.path(data_dir, "pre_5_lf_vr.xlsx")),
  IRL = list(cop = file.path(data_dir, "pre_5_cop_irl.xlsx"),
             lf  = file.path(data_dir, "pre_5_lf_irl.xlsx"))
)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl first: install.packages('readxl')")
}
missing <- unlist(files, use.names = FALSE)
missing <- missing[!file.exists(missing)]
if (length(missing)) stop("Missing workbook(s):\n", paste(missing, collapse = "\n"))
if (any(trial_indices < 1 | trial_indices > 25) || anyDuplicated(trial_indices)) {
  stop("trial_indices must be unique row indices between 1 and 25.")
}

# Read participant IDs from the first row and trials from rest of rows
read_wide <- function(path) {
  workbook <- as.data.frame(readxl::read_excel(path, col_names = FALSE),
                            check.names = FALSE)
  if (nrow(workbook) < 26L) stop("Expected one ID row plus 25 data rows: ", path)
  ids <- trimws(as.character(unlist(workbook[1, , drop = FALSE], use.names = FALSE)))
  keep <- !is.na(ids) & nzchar(ids)
  values <- workbook[2:26, keep, drop = FALSE]
  values[] <- lapply(values, function(x) suppressWarnings(as.numeric(as.character(x))))
  names(values) <- ids[keep]
  if (anyDuplicated(names(values))) stop("Duplicate participant IDs in: ", path)
  values
}

# Calculate an approximate 95% confidence interval for Pearson r by forming
# the interval in Fisher-z units and back-transforming its limits to r
fisher_r_ci <- function(r, n, confidence = .95) {
  if (!is.finite(r) || n <= 3L || abs(r) >= 1) return(c(NA_real_, NA_real_))
  z <- atanh(r)
  se_z <- 1 / sqrt(n - 3)
  critical <- qnorm(1 - (1 - confidence) / 2)
  tanh(z + c(-1, 1) * critical * se_z)
}

analyze_environment <- function(environment, cop_wide, lf_wide) {
  ids <- intersect(names(cop_wide), names(lf_wide))
  if (!length(ids)) stop("No matching COP/LF participant columns for ", environment)
  unmatched <- union(setdiff(names(cop_wide), names(lf_wide)),
                     setdiff(names(lf_wide), names(cop_wide)))
  if (length(unmatched)) {
    warning("Excluding unmatched participant column(s) for ", environment, ": ",
            paste(unmatched, collapse = ", "))
  }
  
  participant_results <- lapply(ids, function(id) {
    subject_number <- suppressWarnings(as.integer(id))
    if (is.na(subject_number)) {
      stop("Participant IDs must be numeric for the even-ID sign-flip rule: ", id)
    }
    
    # Select the same trial rows from COPdiff and LFdiff
    cop <- as.numeric(cop_wide[[id]])[trial_indices]
    lf <- as.numeric(lf_wide[[id]])[trial_indices]
    
    # The study's existing coding treats even-numbered participants as R-to-L starters
    if (subject_number %% 2L == 0L) {
      cop <- -cop
      lf <- -lf
    }
    
    # -999 excluded
    valid <- is.finite(cop) & is.finite(lf) & cop != -999 & lf != -999
    cop <- cop[valid]
    lf <- lf[valid]
    
    r <- p <- z <- NA_real_
    ci <- c(NA_real_, NA_real_)
    n <- length(cop)
    if (n >= 3L && sd(cop) > 0 && sd(lf) > 0) {
      correlation_test <- cor.test(cop, lf, method = "pearson", conf.level = .95)
      r <- unname(correlation_test$estimate)
      p <- correlation_test$p.value
      z <- atanh(r) # Fisher r-to-z: z = atanh(r)
      ci <- fisher_r_ci(r, n) # 95% CI for r using Fisher z, transformed back
    }
    
    list(id = id, subject_number = subject_number, n = n,
         r = r, ci_low = ci[1], ci_high = ci[2], z = z, p = p)
  })
  
  list(environment = environment, participants = participant_results)
}

cop <- lapply(files, function(x) read_wide(x$cop))
lf <- lapply(files, function(x) read_wide(x$lf))
results <- lapply(names(files), function(env) {
  analyze_environment(env, cop[[env]], lf[[env]])
})
names(results) <- names(files)

# Assemble a participant-by-environment table
participant_table <- do.call(rbind, lapply(results, function(result) {
  do.call(rbind, lapply(result$participants, function(x) {
    data.frame(Subject = x$id, Environment = result$environment,
               n = x$n, r = x$r, CI95_low = x$ci_low,
               CI95_high = x$ci_high, Fisher_z = x$z, p = x$p,
               stringsAsFactors = FALSE)
  }))
}))
wide_z <- reshape(participant_table[, c("Subject", "Environment", "Fisher_z")],
                  idvar = "Subject", timevar = "Environment", direction = "wide")
paired <- wide_z[is.finite(wide_z$Fisher_z.VR) & is.finite(wide_z$Fisher_z.IRL), ]
if (nrow(paired) < 2L) stop("At least two participants need valid correlations in both environments.")
paired$z_difference_VR_minus_IRL <- paired$Fisher_z.VR - paired$Fisher_z.IRL
paired_test <- t.test(paired$Fisher_z.VR, paired$Fisher_z.IRL,
                      paired = TRUE, conf.level = .95)

# Print participant-level correlations, each Fisher transformation, and the
# back-transformation
# Group mean r is tanh(mean participant Fisher z)
cat("Analysis window:", window_name, "\n")
cat("Measurement rows used:", paste(trial_indices, collapse = ", "), "\n")
cat("Missing-value code excluded: -999\n\n")
cat("Participant-level correlations and Fisher transformations\n")
cat("For each participant: z = atanh(r); back-transformed r = tanh(z).\n")
print(participant_table, row.names = FALSE, digits = 4)

cat("\nEnvironment-level summaries\n")
for (env in names(results)) {
  zvals <- vapply(results[[env]]$participants, function(x) x$z, numeric(1))
  zvals <- zvals[is.finite(zvals)]
  if (length(zvals) < 2L) {
    cat(env, ": fewer than two valid participant correlations; group CI unavailable.\n")
    next
  }
  mean_z <- mean(zvals)
  mean_z_ci <- t.test(zvals, conf.level = .95)$conf.int
  mean_r <- tanh(mean_z)
  mean_r_ci <- tanh(mean_z_ci)
  cat(sprintf(paste0(
    "%s: mean Fisher z = %.3f, 95%% CI [%.3f, %.3f]; ",
    "back-transformed group mean r = tanh(%.3f) = %.3f, ",
    "95%% CI [%.3f, %.3f]; n participants = %d\n"),
    env, mean_z, mean_z_ci[1], mean_z_ci[2], mean_z, mean_r,
    mean_r_ci[1], mean_r_ci[2], length(zvals)))
}

cat("\nPaired environment comparison of participant Fisher z values\n")
cat(sprintf(paste0("VR minus IRL: mean difference in z = %.3f, ",
                   "95%% CI [%.3f, %.3f], t(%g) = %.3f, p = %s, n = %d paired participants\n"),
            mean(paired$z_difference_VR_minus_IRL),
            paired_test$conf.int[1], paired_test$conf.int[2],
            unname(paired_test$parameter), unname(paired_test$statistic),
            format.pval(paired_test$p.value, digits = 4), nrow(paired)))
