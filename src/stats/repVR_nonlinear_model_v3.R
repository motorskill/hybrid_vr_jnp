###############################################################
# repVR: nonlinear mixed-effects model of pre-switch MCom
# Reanalysis for Reviewer 2 (comments 6 & 7)
#
# HOW TO RUN
#   1. Put this script and the two data files in the SAME folder:
#        pre_five_vr.xlsx
#        pre_5_irl.xlsx
#   2. Open this script in RStudio.
#   3. Menu: Session > Set Working Directory > To Source File Location
#   4. First time only, install packages (remove the # and run the line):
#        # install.packages(c("openxlsx", "nlme", "emmeans"))
#   5. Click "Source" (or select all and Run). It takes a few minutes.
#      Each block of output is headed with ===== and says which
#      manuscript table/row it supplies.
###############################################################

library(openxlsx)
library(nlme)
library(emmeans)

###############################################################
# 1. LOAD AND RESHAPE THE DATA
###############################################################
read_env <- function(file, env) {
  x <- read.xlsx(file, colNames = TRUE)
  x <- x[, colSums(is.na(x)) < nrow(x), drop = FALSE]      # drop empty columns
  # one column per subject, 25 rows (5 blocks x 5 trials) -> long format
  do.call(rbind, lapply(names(x), function(s)
    data.frame(GlobalTrial = seq_len(nrow(x)), Subject = s,
               MCom = x[[s]], Environment = env)))
}
long <- rbind(read_env("pre_five_vr.xlsx", "VR"),
              read_env("pre_5_irl.xlsx",   "IRL"))

long$Trial   <- ((long$GlobalTrial - 1) %% 5) + 1           # trial within block, 1-5
long$Block   <- ((long$GlobalTrial - 1) %/% 5) + 1          # block 1-5
long$Trial_c <- long$Trial - 1                              # t in the model, 0-4

VR_first <- c("1", "2", "3", "4", "9", "11", "13", "15")    # participants who did VR first
long$Session <- ifelse((long$Subject %in% VR_first) == (long$Environment == "VR"),
                       "Session1", "Session2")
long$Order   <- ifelse(long$Subject %in% VR_first, "VRfirst", "IRLfirst")

long$Subject     <- factor(long$Subject)
long$Environment <- factor(long$Environment, levels = c("IRL", "VR"))
long$Session     <- factor(long$Session, levels = c("Session1", "Session2"))
long$Block       <- factor(long$Block)

# 999 = missing-data code (Subject 4, real-world, Block 4, Trial 2; Excel row 18 of pre_5_irl.xlsx).
# Treat 999 and any blank cell as missing, then drop those trials before fitting.
long$MCom[long$MCom == 999] <- NA
cat("\n===== Trials treated as missing (999 or blank) =====\n")
print(long[is.na(long$MCom), c("Subject", "Environment", "Block", "Trial", "MCom")])
long <- long[!is.na(long$MCom), ]
cat("Observations:", nrow(long), "   Participants:", nlevels(long$Subject),
    "   Max MCom used:", max(long$MCom), "\n")
stopifnot(nrow(long) == 749, !any(long$MCom == 999))   # stop if the missing trial was not removed

###############################################################
# 2. SETTINGS
###############################################################
# Sum-to-zero (effect) coding: main effects = averaged over the other factors
options(contrasts = c("contr.sum", "contr.poly"))
emm_options(lmer.df = "asymptotic")

ctl <- nlmeControl(maxIter = 200, pnlsMaxIter = 50, msMaxIter = 200)

# Model:  MCom = A - B * exp(-k * t)
#   A = asymptote (Nmm), B = beta, improvement from trial 1 to asymptote (Nmm), k = learning rate

# Number of k coefficients for the primary model (1 intercept + 17 effects = 18)
nk <- ncol(model.matrix(~ Block * Session + Environment * Session + Block:Environment, long))

###############################################################
# 3. PRIMARY MODEL
#    k ~ Block*Session + Environment*Session + Block:Environment
###############################################################
# Primary model: random effect on A (participant-specific asymptote).
# A random effect on k was also tested (random = pdDiag(A + k ~ 1)). With the corrected data this
# model did not converge (also not with 2000 iterations, alternative optimizers, or k on the log
# scale), so it is not used for inference. Its estimates at termination gave the same conclusions
# (no significant effects on k; all p >= .09) -- see Methods and Response 6.2.
fit <- nlme(
  MCom ~ A - B * exp(-k * Trial_c), data = long,
  fixed  = list(A ~ 1, B ~ 1,
                k ~ Block * Session + Environment * Session + Block:Environment),
  random = pdDiag(A ~ 1), groups = ~ Subject,
  start  = c(mean(long$MCom), 100, 0.5, rep(0, nk - 1)),
  control = ctl)
fit_startvals <- fit   # same model; used below as the comparison model for the asymptote LRT

cat("\n===== Random-effect SDs  -> Results (last sentence), Table S1 'Random effects', Response 6.2 =====\n")
print(VarCorr(fit))       # StdDev column: A (Nmm), Residual (Nmm)

cat("\n===== A, beta, mean k with 95% CI  -> Results (first sentence), Table S1 'Parameters' =====\n")
print(intervals(fit, which = "fixed")$fixed[c("A", "B", "k.(Intercept)"), ])

cat("\n===== F-tests on k  -> Results, Table S1 'Effects on learning rate', Response 6.3 (report rows starting k.) =====\n")
print(anova(fit, type = "marginal"))

cat("\n===== Environment: hybrid-VR - real-world dk  -> Results, Table S1, Response 6.3 =====\n")
print(confint(pairs(emmeans(fit, ~ Environment, param = "k"), reverse = TRUE)))

cat("\n===== Session: Session 2 - Session 1 dk  -> Results, Table S1, Response 6.3 =====\n")
print(confint(pairs(emmeans(fit, ~ Session, param = "k"), reverse = TRUE)))

cat("\n===== Environment x Session dk  -> Results, Table S1, Response 6.3 =====\n")
print(confint(contrast(emmeans(fit, ~ Environment * Session, param = "k"),
                       interaction = c("revpairwise", "revpairwise"))))

cat("\n===== Supplementary Table 2, Real-world / Hybrid-VR columns: k by Block x Environment =====\n")
print(emmeans(fit, ~ Environment | Block, param = "k"))   # emmean, lower.CL, upper.CL

cat("\n===== Supplementary Table 2, Session 1 / Session 2 columns: k by Block x Session =====\n")
print(emmeans(fit, ~ Session | Block, param = "k"))       # emmean, lower.CL, upper.CL

###############################################################
# 4. ASYMPTOTE (A) CHECKS: does level differ by Environment / Session / Block?
#    Random effect on A only, as in the primary model
#    Likelihood-ratio test:  chi2 = 2 x (logLik bigger - logLik smaller)
###############################################################
lrt <- function(smaller, bigger) {
  chi2 <- 2 * as.numeric(logLik(bigger) - logLik(smaller))
  df   <- attr(logLik(bigger), "df") - attr(logLik(smaller), "df")
  c(chi2 = round(chi2, 2), df = df, p = round(pchisq(chi2, df, lower.tail = FALSE), 3))
}

# A ~ Environment * Session  (A gets 4 coefficients)
fit_A1 <- nlme(
  MCom ~ A - B * exp(-k * Trial_c), data = long,
  fixed  = list(A ~ Environment * Session, B ~ 1,
                k ~ Block * Session + Environment * Session + Block:Environment),
  random = pdDiag(A ~ 1), groups = ~ Subject,
  start  = c(mean(long$MCom), 0, 0, 0, 100, 0.5, rep(0, nk - 1)),
  control = ctl)

cat("\n===== LRT: Environment, Session, Env x Session on A  -> Results, Table S1 'Asymptote', Response 6.4 =====\n")
print(lrt(fit_startvals, fit_A1))

cat("\n===== Env x Session (order) contrast on A, df = 13  -> Table S1 'Asymptote', Response 6.4 =====\n")
emA <- emmeans(fit_A1, ~ Environment * Session, param = "A")
print(confint(contrast(emA, interaction = c("revpairwise", "revpairwise")), df = 13))
print(test(contrast(emA, interaction = c("revpairwise", "revpairwise")), df = 13))

# Does the Env x Session difference in A vary across blocks? (Response 6.4)
# Higher-order A terms do not converge with Block:Environment on k,
# so they are tested with the reduced k model  k ~ Block*Session + Environment*Session
nk0 <- ncol(model.matrix(~ Block * Session + Environment * Session, long))
fit_R3 <- nlme(MCom ~ A - B * exp(-k * Trial_c), data = long,
  fixed = list(A ~ (Environment + Session) * Block + Environment:Session, B ~ 1, k ~ Block * Session + Environment * Session),
  random = pdDiag(A ~ 1), groups = ~ Subject,
  start = c(mean(long$MCom), rep(0, 15), 100, 0.5, rep(0, nk0 - 1)), control = ctl)
fit_R4 <- nlme(MCom ~ A - B * exp(-k * Trial_c), data = long,
  fixed = list(A ~ Environment * Session * Block, B ~ 1, k ~ Block * Session + Environment * Session),
  random = pdDiag(A ~ 1), groups = ~ Subject,
  start = c(mean(long$MCom), rep(0, 19), 100, 0.5, rep(0, nk0 - 1)), control = ctl)

cat("\n===== LRT: does the Env x Session difference in A vary across blocks?  -> Response 6.4 =====\n")
print(lrt(fit_R3, fit_R4))

###############################################################
# 5. DESCRIPTIVES FOR RESPONSE 6.4: mean MCom on trials 3-5 by order group
###############################################################
cat("\n===== Group means (Nmm), trials 3-5  -> Response 6.4 =====\n")
subj <- aggregate(MCom ~ Subject + Order + Environment + Session, long[long$Trial >= 3, ], mean)
print(aggregate(MCom ~ Order + Environment + Session, subj, function(x) round(mean(x), 0)))

###############################################################
# 6. FULL THREE-WAY ON k: does not converge  -> Methods, Response 6.3
###############################################################
nk3 <- ncol(model.matrix(~ Block * Session * Environment, long))
fit_3way <- try(nlme(MCom ~ A - B * exp(-k * Trial_c), data = long,
  fixed = list(A ~ 1, B ~ 1, k ~ Block * Session * Environment),
  random = pdDiag(A ~ 1), groups = ~ Subject,
  start = c(mean(long$MCom), 100, 0.5, rep(0, nk3 - 1)), control = ctl))
cat("\n===== Three-way model on k (expected to fail) =====\n")
print(class(fit_3way))

###############################################################
# 7. MANUSCRIPT-READY NUMBERS
#    Every reported value is printed here already formatted.
#    Copy from this output (or report_values.csv) -- never retype or round by hand.
###############################################################
# --- formatting helpers ---
fmt_num <- function(x, d = 2) {                    # fixed decimals, true minus sign
  s <- formatC(x, format = "f", digits = d)
  s <- sub("^-", "\u2212", s)
  sub("^\u2212(0\\.?0*)$", "\\1", s)       # no "-0.00"
}
fmt_p  <- function(p) ifelse(p < .001, "p < .001", paste0("p = ", sub("^0", "", formatC(p, format = "f", digits = 3))))
fmt_ci <- function(est, lo, hi, d = 2) sprintf("%s [%s, %s]", fmt_num(est, d), fmt_num(lo, d), fmt_num(hi, d))
fmt_F  <- function(tab, row) sprintf("F(%d, %d) = %s, %s", tab[row, "numDF"], tab[row, "denDF"],
                                     fmt_num(tab[row, "F-value"]), fmt_p(tab[row, "p-value"]))
fmt_contrast <- function(ct, d = 2) { s <- as.data.frame(summary(ct, infer = TRUE))
  sprintf("%s, 95%% CI [%s, %s]", fmt_num(s$estimate, d), fmt_num(s$lower.CL, d), fmt_num(s$upper.CL, d)) }

# --- collect every reported value in one table ---
aov_k <- anova(fit, type = "marginal")
ivals <- intervals(fit, which = "fixed")$fixed
vals <- c(
  "beta (Nmm)"                    = fmt_ci(ivals["B", 2], ivals["B", 1], ivals["B", 3], 0),
  "A asymptote (Nmm)"             = fmt_ci(ivals["A", 2], ivals["A", 1], ivals["A", 3], 0),
  "k mean"                        = fmt_ci(ivals["k.(Intercept)", 2], ivals["k.(Intercept)", 1], ivals["k.(Intercept)", 3]),
  "Environment F"                 = fmt_F(aov_k, "k.Environment"),
  "Environment dk"                = fmt_contrast(pairs(emmeans(fit, ~ Environment, param = "k"), reverse = TRUE)),
  "Session F"                     = fmt_F(aov_k, "k.Session"),
  "Session dk"                    = fmt_contrast(pairs(emmeans(fit, ~ Session, param = "k"), reverse = TRUE)),
  "Environment x Session F"       = fmt_F(aov_k, "k.Session:Environment"),
  "Environment x Session dk"      = fmt_contrast(contrast(emmeans(fit, ~ Environment * Session, param = "k"),
                                                          interaction = c("revpairwise", "revpairwise"))),
  "Block F"                       = fmt_F(aov_k, "k.Block"),
  "Block x Session F"             = fmt_F(aov_k, "k.Block:Session"),
  "Block x Environment F"         = fmt_F(aov_k, "k.Block:Environment"),
  "SD A (Nmm)"                    = fmt_num(as.numeric(VarCorr(fit)["A", "StdDev"]), 0),
  "SD residual (Nmm)"             = fmt_num(fit$sigma, 0)
)
lr1 <- lrt(fit_startvals, fit_A1); lr2 <- lrt(fit_R3, fit_R4)
vals["Asymptote LRT"]          <- sprintf("\u03c7\u00b2(%d) = %s, %s", lr1[["df"]], fmt_num(lr1[["chi2"]]), fmt_p(lr1[["p"]]))
vals["Asymptote x Block LRT"]  <- sprintf("\u03c7\u00b2(%d) = %s, %s", lr2[["df"]], fmt_num(lr2[["chi2"]]), fmt_p(lr2[["p"]]))
cA <- as.data.frame(summary(contrast(emA, interaction = c("revpairwise", "revpairwise")), df = 13, infer = TRUE))
vals["Asymptote Env x Session"] <- sprintf("%s Nmm, 95%% CI [%s, %s], %s", fmt_num(cA$estimate, 0), fmt_num(cA$lower.CL, 0),
                                          fmt_num(cA$upper.CL, 0), fmt_p(cA$p.value))

# Supplementary Table 2 cells
s2e <- as.data.frame(summary(emmeans(fit, ~ Environment | Block, param = "k")))
s2s <- as.data.frame(summary(emmeans(fit, ~ Session | Block, param = "k")))
for (i in seq_len(nrow(s2e))) vals[sprintf("Table S2 Block %s %s", s2e$Block[i], s2e$Environment[i])] <- fmt_ci(s2e$emmean[i], s2e$lower.CL[i], s2e$upper.CL[i])
for (i in seq_len(nrow(s2s))) vals[sprintf("Table S2 Block %s %s", s2s$Block[i], s2s$Session[i])]     <- fmt_ci(s2s$emmean[i], s2s$lower.CL[i], s2s$upper.CL[i])

cat("\n===== MANUSCRIPT-READY VALUES (copy from here or from report_values.csv) =====\n")
out <- data.frame(item = names(vals), value = unname(vals))
cat(sprintf("%-28s %s\n", out$item, out$value), sep = "")
write.csv(out, "report_values.csv", row.names = FALSE, fileEncoding = "UTF-8")
