###############################################################################
# Integrated mortality analyses for:
# Bumblebee_pesticide_susceptibility_and_rearing_environment.csv
#
# Purpose
# -------
# This script reproduces the analyses and figures for Experiments 1–4 in the
# linked manuscript script examining pesticide sensitivity in wild and 
# lab reared bumblebees
#
# Global paper-wide inclusion rules
# ---------------------------------
#   - Bombus terrestris only
#   - Worker caste only
#   - Males excluded
#   - Valid mortality data only
#
# General structure within each experiment
# ----------------------------------------
#   A. Filter the consolidated dataset to the manuscript analysis set
#   B. Audit sample sizes and key variables after filtering
#   C. Inspect survival patterns with Kaplan–Meier curves
#   D. Run the manuscript-relevant inferential statistics
#   E. Recreate the publication figure
#   F. Explore supporting explanatory variables (e.g. weight, sucrose)
###############################################################################


###############################################################################
# 0. HOUSEKEEPING AND DATA IMPORT
###############################################################################

# Start from a clean R session and close open graphics devices so that previous
# objects or figures do not interfere with this run.
rm(list = ls())
graphics.off()

# Load packages.
# survival / survminer / coxme are used for survival analyses.
# ggplot2 is used for custom figures.
# car, e1071, FSA are used for diagnostics and non-parametric follow-up tests.
# MuMIn and emmeans are used for model comparison and post hoc contrasts.
library(survival)
library(ggplot2)
library(survminer)
library(car)
library(e1071)
library(FSA)
library(coxme)
library(MuMIn)
library(emmeans)
library(ggfortify)

# Load the consolidated dataset.
Guide <- read.csv(
  file = "Bumblebee_pesticide_susceptibility_and_rearing_environment.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Quick visual check that the dataset imported correctly.
head(Guide)
tail(Guide)

# Print column names to confirm the consolidated dataset structure.
print(colnames(Guide))

# Check that the minimum required columns for the integrated workflow are present.
required_cols <- c(
  "Experiment", "Experimental_batch", "Caste", "Treatment", "Species", "Male",
  "Bee_weight_g", "Time_in_experiment_hours", "Mortality", "Valid_mort_data"
)
missing_cols <- setdiff(required_cols, names(Guide))
if (length(missing_cols) > 0) {
  stop(paste(
    "The following required columns are missing from the consolidated dataset:",
    paste(missing_cols, collapse = ", ")
  ))
}

# Standardise key variable types. This prevents downstream problems caused by
# inconsistent import of numeric or categorical columns as character/factor.
Guide$Experiment <- trimws(as.character(Guide$Experiment))
Guide$Caste <- trimws(as.character(Guide$Caste))
Guide$Treatment <- trimws(as.character(Guide$Treatment))
Guide$Species <- trimws(as.character(Guide$Species))
Guide$Male <- suppressWarnings(as.numeric(as.character(Guide$Male)))
Guide$Valid_mort_data <- suppressWarnings(as.numeric(as.character(Guide$Valid_mort_data)))
Guide$Bee_weight_g <- suppressWarnings(as.numeric(as.character(Guide$Bee_weight_g)))
Guide$Time_in_experiment_hours <- suppressWarnings(as.numeric(as.character(Guide$Time_in_experiment_hours)))
Guide$Mortality <- suppressWarnings(as.numeric(as.character(Guide$Mortality)))


###############################################################################
# EXPERIMENT 1
# Wild-caught vs lab-reared commercial workers
# Focal pesticide: thiamethoxam
# Positive control: dimethoate
#
# Manuscript logic
# ----------------
#   - Survival curves are used to visualise mortality over time.
#   - Positive control is retained for plotting and context.
#   - Positive control is excluded from inferential mortality comparisons.
#   - A Cox model is attempted first, but the manuscript analysis uses
#     chi-squared / Fisher tests because mortality is sparse or highly separated.
###############################################################################

###############################################################################
# 1A. Filter the consolidated dataset to the Experiment 1 analysis set
###############################################################################

# Explicitly apply the paper-wide inclusion rules to Experiment 1.
Guide1 <- subset(
  Guide,
  Experiment == "Experiment01" &
    Caste == "Worker" &
    Species == "terrestris" &
    Male != 1 &
    Valid_mort_data == 1
)

# Force treatment order so that legends and survival curves follow the intended
# manuscript figure layout.
Guide1$Treatment <- factor(Guide1$Treatment, levels = c("RC", "WC", "WP", "RP", "DIM"))

# Colour palette used throughout Experiment 1 figure production.
trt_cols <- c(
  "RC"  = "#56B4E9",  # Reared Control
  "WC"  = "#F0E442",  # Wild Control
  "WP"  = "#009E73",  # Wild Pesticide
  "RP"  = "#E69F00",  # Reared Pesticide
  "DIM" = "#000000"   # Positive Control
)

###############################################################################
# 1B. Audit sample sizes and key variables after filtering
###############################################################################

# These checks confirm that the Experiment 1 subset contains the expected rows
# and that the global inclusion criteria have been applied correctly.
exp01_count <- sum(Guide$Experiment == "Experiment01", na.rm = TRUE)
print(exp01_count)

worker_count_exp01 <- sum(Guide$Experiment == "Experiment01" & Guide$Caste == "Worker", na.rm = TRUE)
print(worker_count_exp01)

lucorum_count_exp01 <- sum(Guide$Experiment == "Experiment01" & Guide$Species == "lucorum", na.rm = TRUE)
print(lucorum_count_exp01)

terrestris_count_exp01 <- sum(Guide$Experiment == "Experiment01" & Guide$Species == "terrestris", na.rm = TRUE)
print(terrestris_count_exp01)

male_count_exp01 <- sum(Guide$Experiment == "Experiment01" & Guide$Male == 1, na.rm = TRUE)
print(male_count_exp01)

valid_mort_count_exp01 <- sum(Guide$Experiment == "Experiment01" & Guide$Valid_mort_data == 1, na.rm = TRUE)
print(valid_mort_count_exp01)

# Final treatment counts after all filtering.
print(table(Guide1$Treatment, useNA = "ifany"))

# Expected manuscript sample sizes:
#   RC = 45, RP = 39, WC = 46, WP = 45, DIM = 19

# Quick audit of core variables needed for modelling and plotting.
sum(is.na(Guide1$Time_in_experiment_hours))
sum(is.na(Guide1$Mortality))
sum(is.na(Guide1$Bee_weight_g))

###############################################################################
# 1C. Survival inspection: Kaplan–Meier curves
###############################################################################

# Fit survival curves to inspect whether overall mortality patterns match the
# original analysis script and the manuscript figure.
Guide1$Treatment <- factor(Guide1$Treatment, levels = c("RC", "WC", "WP", "RP", "DIM"))
km_trt_fit <- survfit(Surv(Time_in_experiment_hours, Mortality) ~ Treatment, data = Guide1)
summary(km_trt_fit)

# Quick diagnostic plots used during development to verify:
#   (i) treatment coding is correct,
#   (ii) all expected groups are present,
#   (iii) survival patterns match the original script.
ggsurvplot(km_trt_fit)
ggsurvplot(km_trt_fit, pval = TRUE)
ggsurvplot(km_trt_fit, pval = FALSE, legend = "bottom")

###############################################################################
# 1D. Mortality inference for the manuscript
###############################################################################

# Positive control is excluded from formal treatment-vs-treatment inference.
NoDimGuide <- subset(Guide1, Treatment != "DIM")
NoDimGuide$Treatment <- droplevels(NoDimGuide$Treatment)

# A Cox proportional hazards model is attempted first because it preserves the
# time-to-event information, but Experiment 1 mortality patterns are too
# separated/sparse for the fitted model to be useful.
PH <- coxph(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    as.factor(Treatment) + Bee_weight_g,
  data = NoDimGuide
)
summary(PH)

# The manuscript-facing analysis therefore uses contingency-table methods.
# Chi-squared is used as the default global and pairwise test.
# Fisher’s exact test is used as a robustness check when expected counts are low.
NoDIMtbl <- table(NoDimGuide$Treatment, NoDimGuide$Mortality)
NoDIMtbl
# 0 = survived / censored at 48 h
# 1 = died

NoDIMcs <- chisq.test(NoDIMtbl)
NoDIMcs
# Global treatment effect is expected to be highly significant.

# There are six pairwise comparisons in Experiment 1.
# Bonferroni-adjusted alpha: 0.05 / 6 = 0.008333.

# Helper function:
# After subsetting a factor, unused levels can remain attached. droplevels()
# removes those unused levels so the resulting contingency table contains only
# the two groups being compared.
make_pair_table <- function(data, grp1, grp2) {
  subdat <- droplevels(subset(data, Treatment %in% c(grp1, grp2)))
  table(subdat$Treatment, subdat$Mortality)
}

# 1. RC vs RP
RCRPtbl <- make_pair_table(Guide1, "RC", "RP")
RCRPtbl
RCRPcs <- chisq.test(RCRPtbl)
RCRPcs

# 2. RC vs WC
RCWCtbl <- make_pair_table(Guide1, "RC", "WC")
RCWCtbl
RCWCcs <- chisq.test(RCWCtbl)
RCWCcs
# The counts are very similar here, so Fisher’s exact test is also reported.
RCWCfisher_result <- fisher.test(RCWCtbl)
print(RCWCfisher_result)

# 3. RC vs WP
RCWPtbl <- make_pair_table(Guide1, "RC", "WP")
RCWPtbl
RCWPcs <- chisq.test(RCWPtbl)
RCWPcs
# Again, this comparison has sparse counts, so Fisher’s exact test is used as a
# robustness check.
RCWPfisher_result <- fisher.test(RCWPtbl)
print(RCWPfisher_result)

# 4. RP vs WC
RPWCtbl <- make_pair_table(Guide1, "RP", "WC")
RPWCtbl
RPWCcs <- chisq.test(RPWCtbl)
RPWCcs

# 5. RP vs WP
RPWPtbl <- make_pair_table(Guide1, "RP", "WP")
RPWPtbl
RPWPcs <- chisq.test(RPWPtbl)
RPWPcs

# 6. WC vs WP
WCWPtbl <- make_pair_table(Guide1, "WC", "WP")
WCWPtbl
WCWPcs <- chisq.test(WCWPtbl)
WCWPcs
# Sparse counts again motivate a Fisher check.
WCWPfisher_result <- fisher.test(WCWPtbl)
print(WCWPfisher_result)

# Summary of manuscript-relevant pairwise outputs.
RCRPcs
RCWCcs
print(RCWCfisher_result)
RCWPcs
print(RCWPfisher_result)
RPWCcs
RPWPcs
WCWPcs
print(WCWPfisher_result)

###############################################################################
# 1E. Figure 1
#
# Purpose
# -------
# Recreate the Experiment 1 manuscript survival figure.
#
# Notes
# -----
#   - Survival is plotted on the native 0–1 scale and displayed as percentages
#     using surv.scale = "percent".
#   - Significance letters are added manually to match the manuscript layout.
#   - Letter positions are hand-tuned and may need small adjustments if export
#     dimensions are changed.
###############################################################################

Guide1$Treatment <- factor(Guide1$Treatment, levels = c("RC", "WC", "WP", "RP", "DIM"))
km_trt_fit <- survfit(Surv(Time_in_experiment_hours, Mortality) ~ Treatment, data = Guide1)
summary(km_trt_fit)

# Quick visual checks used during figure construction.
ggsurvplot(km_trt_fit)
ggsurvplot(km_trt_fit, pval = TRUE)
ggsurvplot(km_trt_fit, pval = FALSE, legend = "bottom")

a <- ggsurvplot(
  km_trt_fit,
  data = Guide1,
  size = 1.8,
  legend = "right",
  legend.labs = c(
    "Reared Control",
    "Wild Control",
    "Wild Pesticide",
    "Reared Pesticide",
    "Positive Control"
  ),
  legend.title = "",
  font.legend = 18,
  font.x = 26,
  font.y = 26,
  font.xtickslab = 18,
  font.ytickslab = 18,
  break.time.by = 24,
  xlim = c(0, 50),
  ylim = c(0, 1.05),
  censor = FALSE,
  conf.int = FALSE,
  xlab = "Time (hours)",
  ylab = "Survival",
  surv.scale = "percent",
  palette = unname(trt_cols[c("RC", "WC", "WP", "RP", "DIM")]),
  ggtheme = theme_classic(base_size = 18) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(colour = "black"),
      legend.position = "right",
      legend.text = element_text(size = 18, colour = "black"),
      legend.key = element_blank(),
      plot.margin = margin(t = 10, r = 65, b = 10, l = 10)
    )
)

# The manuscript figure uses significance letters at 48 h:
#   RC, WC and WP = a
#   RP            = b
#   DIM           = c
km_48 <- summary(km_trt_fit, times = 48, extend = TRUE)

letters_df <- data.frame(
  Treatment = gsub("^Treatment=", "", km_48$strata),
  survival = km_48$surv,
  stringsAsFactors = FALSE
)

letters_df <- letters_df[match(c("RC", "WC", "WP", "RP", "DIM"), letters_df$Treatment), ]
letters_df$label <- c("a", "a", "a", "b", "c")

# Annotation coordinates are on the raw 0–1 scale because surv.scale = "percent"
# only affects display.
letters_df$x <- c(49.3, 49.3, 49.3, 49.3, 5.2)
letters_df$y <- c(
  letters_df$survival[letters_df$Treatment == "RC"] + 0.005,
  letters_df$survival[letters_df$Treatment == "WC"] - 0.001,
  letters_df$survival[letters_df$Treatment == "WP"] - 0.002,
  letters_df$survival[letters_df$Treatment == "RP"] + 0.007,
  letters_df$survival[letters_df$Treatment == "DIM"] + 0.015
)

letters_df$colour <- trt_cols[letters_df$Treatment]

for (i in seq_len(nrow(letters_df))) {
  a$plot <- a$plot +
    annotate(
      "text",
      x = letters_df$x[i],
      y = letters_df$y[i],
      label = letters_df$label[i],
      colour = letters_df$colour[i],
      size = 8,
      hjust = 0,
      family = "sans"
    )
}

print(a$plot)

ggsave(
  filename = "Experiment01_Fig1_wild_vs_reared_bt_worker_mortality.pdf",
  plot = a$plot,
  width = 11,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Experiment01_Fig1_wild_vs_reared_bt_worker_mortality.png",
  plot = a$plot,
  width = 11,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

###############################################################################
# 1F. Supporting explanatory variables: weight and sucrose
###############################################################################

# Check whether body weight differs by mortality status or treatment.
plot(Bee_weight_g ~ Mortality, data = NoDimGuide)
plot(Bee_weight_g ~ as.factor(Treatment), data = NoDimGuide)

# Treatment-level differences in body weight.
WeightComp <- aov(Bee_weight_g ~ as.factor(Treatment), data = NoDimGuide)
summary(WeightComp)
WeightCompTukey_Results <- TukeyHSD(WeightComp)
print(WeightCompTukey_Results)
# Expected interpretation:
# both wild treatments are heavier than both reared treatments.

# Mortality is concentrated in RP, so within-treatment weight effects are only
# meaningfully testable there.
RPOnlyGuide <- droplevels(subset(NoDimGuide, Treatment == "RP"))

plot(Bee_weight_g ~ Mortality, data = RPOnlyGuide)
plot(Bee_weight_g ~ as.factor(Treatment), data = RPOnlyGuide)

# Preserve the original modelling direction used in the legacy script:
# weight is the response, mortality class is the predictor.
RPMortVSWeight <- glm(Bee_weight_g ~ as.factor(Mortality), data = RPOnlyGuide)
summary(RPMortVSWeight)
confint(RPMortVSWeight)

# Examine the RP vs WP comparison directly.
RPWPGuide <- droplevels(subset(Guide1, Treatment %in% c("RP", "WP")))
head(RPWPGuide)

f <- ggplot(
  data = RPWPGuide,
  aes(y = as.factor(Mortality), x = Bee_weight_g, color = Treatment)
) +
  geom_point(position = position_jitter(width = 0, height = 0.1), size = 2) +
  geom_vline(xintercept = 0.145, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.204, linetype = "dashed", color = "red") +
  theme_bw() +
  coord_flip() +
  scale_x_continuous(limits = c(0, NA)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20)
  ) +
  labs(x = "Bee Weight at Start (g)", y = "Survival") +
  scale_color_manual(values = c("WP" = "#56B4E9", "RP" = "#000000"))
f

g <- f + scale_y_discrete(labels = c("0" = "Survived", "1" = "Died"))
g

ggsave("wild vs reared bt worker weight RPWP scatter.pdf", g, width = 6, height = 8, dpi = 300)

# Match the weight ranges of dead RP bees and alive WP bees to test whether the
# treatment effect persists within a comparable size window.
deadRP <- subset(RPWPGuide, Treatment == "RP" & Mortality == 1)
max(deadRP$Bee_weight_g, na.rm = TRUE)
min(deadRP$Bee_weight_g, na.rm = TRUE)

aliveWP <- subset(RPWPGuide, Treatment == "WP" & Mortality == 0)
max(aliveWP$Bee_weight_g, na.rm = TRUE)
min(aliveWP$Bee_weight_g, na.rm = TRUE)

overlapWPRP <- droplevels(
  subset(
    Guide1,
    (Treatment == "RP" & Bee_weight_g >= 0.145 & Bee_weight_g <= 0.204) |
      (Treatment == "WP" & Bee_weight_g >= 0.145 & Bee_weight_g <= 0.204)
  )
)

max(overlapWPRP$Bee_weight_g, na.rm = TRUE)
min(overlapWPRP$Bee_weight_g, na.rm = TRUE)

f <- ggplot(
  overlapWPRP,
  aes(x = Bee_weight_g, y = as.factor(Mortality), color = Treatment)
) +
  geom_point() +
  theme_bw() +
  coord_flip() +
  scale_x_continuous(limits = c(0, NA)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )
f

overlapRPWPtbl <- table(overlapWPRP$Treatment, overlapWPRP$Mortality)
overlapRPWPtbl

overlapcs <- chisq.test(overlapRPWPtbl)
overlapcs
# Interpretation:
# the reared-vs-wild difference remains significant even when matched on weight.

###############################################################################
# 1G. Sucrose consumption
###############################################################################

SucOnlyGuide <- droplevels(subset(Guide1, Valid_sucrose_data == 1))

# Check that valid sucrose-data filtering behaved as expected.
plot(Valid_sucrose_data ~ as.factor(Treatment), data = SucOnlyGuide)
plot(Sucrose_consumed_per_hour ~ as.factor(Treatment), data = SucOnlyGuide)

# Distributional and model-assumption checks.
hist(
  SucOnlyGuide$Sucrose_consumed_per_hour,
  probability = TRUE,
  main = "Histogram of normal data"
)

shapiro.test(SucOnlyGuide$Sucrose_consumed_per_hour)
bartlett.test(Sucrose_consumed_per_hour ~ Treatment, data = SucOnlyGuide)

mod1 <- lm(Sucrose_consumed_per_hour ~ Treatment, data = SucOnlyGuide)
mean(mod1$residuals)

par(mfrow = c(2, 2))
plot(mod1)
par(mfrow = c(1, 1))

acf(mod1$residuals)
durbinWatsonTest(mod1)
skewness(SucOnlyGuide$Sucrose_consumed_per_hour)

# Due to strong non-normality / skew, use a robust non-parametric approach.
kruskal.test(Sucrose_consumed_per_hour ~ as.factor(Treatment), data = SucOnlyGuide)

PT <- dunnTest(
  Sucrose_consumed_per_hour ~ as.factor(Treatment),
  data = SucOnlyGuide,
  method = "bonferroni"
)
PT

# Check treatment means and standard deviations against source summaries.
mean_values <- aggregate(Sucrose_consumed_per_hour ~ Treatment, data = SucOnlyGuide, mean)
mean_values

sd_values <- aggregate(Sucrose_consumed_per_hour ~ Treatment, data = SucOnlyGuide, sd)
sd_values

treatment_label_map <- c(
  "RC"  = "Reared-Control",
  "RP"  = "Reared-Pesticide",
  "WC"  = "Wild-Control",
  "WP"  = "Wild-Pesticide",
  "DIM" = "Positive-Control"
)

cbPaletteOkabeIto <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#000000")

k <- ggplot(
  data = SucOnlyGuide,
  aes(y = as.factor(Treatment), x = Sucrose_consumed_per_hour, fill = as.factor(Treatment))
) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0, height = 0.1), size = 0.75) +
  theme_bw() +
  scale_fill_manual(values = cbPaletteOkabeIto) +
  coord_flip() +
  scale_x_continuous(limits = c(0, NA)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank()
  ) +
  labs(x = "Sucrose Consumed per Hour Alive (g)", y = "Treatment") +
  labs(fill = "Treatment")
k

l <- k + scale_y_discrete(labels = treatment_label_map)
l

ggsave("wild vs reared bt worker sucrose.pdf", l, width = 6, height = 8, dpi = 300)


###############################################################################
# EXPERIMENT 2
# Wild-caught vs. lab-reared commercial worker bumble bees
# Focal pesticides: Amistar and Cyantraniliprole
# Positive control: Dimethoate
#
# Manuscript logic:
#   1. Build the Experiment 2 analysis subset from the consolidated dataset
#   2. Inspect survival curves for all filtered treatments
#   3. Fit mixed-effects Cox proportional hazards models
#   4. Remove DIM groups for the main manuscript-facing model
#   5. Compare a small set of candidate models and retain the simplest
#      defensible model
#   6. Run six a priori contrasts with Bonferroni correction
#   7. Check proportional hazards and linearity assumptions
#   8. Recreate the publication survival figure (Figure 2)
###############################################################################

###############################################################################
# 2A. Experiment-specific palettes and required columns
###############################################################################

cbPaletteOkabeIto3 <- c(
  "#0072B2",  # Blue
  "#E69F00",  # Orange
  "#56B4E9",  # Sky Blue
  "#009E73",  # Green
  "#F0E442",  # Yellow
  "#CC79A7",  # Pink
  "#000000",  # Black
  "#D55E00",  # Vermilion
  "#999999",  # Gray
  "#9400D3"   # Purple
)

required_cols_exp2 <- c(
  "Experiment", "Experimental_batch", "Caste", "Treatment", "Species", "Male",
  "Bee_weight_g", "Time_in_experiment_hours", "Mortality", "Valid_mort_data",
  "Colony_ID"
)
missing_cols_exp2 <- setdiff(required_cols_exp2, names(Guide))
if (length(missing_cols_exp2) > 0) {
  stop(paste(
    "The following required columns are missing for Experiment 2:",
    paste(missing_cols_exp2, collapse = ", ")
  ))
}

###############################################################################
# 2B. Filter the consolidated dataset to the Experiment 2 analysis set
#
# Paper-wide inclusion rules are applied explicitly:
#   - Experiment02 only
#   - Worker caste only
#   - Bombus terrestris only
#   - Males excluded
#   - Valid mortality data only
###############################################################################

Guide2 <- subset(
  Guide,
  Experiment == "Experiment02" &
    Caste == "Worker" &
    Species == "terrestris" &
    Male != 1 &
    Valid_mort_data == 1
)

# Force treatment order so survival plots, model contrasts, and legends are
# reproducible and aligned with the original Experiment 2 workflow.
Guide2$Treatment <- factor(
  Guide2$Treatment,
  levels = c("ARCON", "RAMI", "RCYN", "WAMI", "WCYN", "WCON", "RDIM", "WDIM")
)

###############################################################################
# 2C. Audit sample sizes and key variables after filtering
###############################################################################

exp02_count <- sum(Guide$Experiment == "Experiment02", na.rm = TRUE)
print(exp02_count)

worker_count_exp02 <- sum(Guide$Experiment == "Experiment02" & Guide$Caste == "Worker", na.rm = TRUE)
print(worker_count_exp02)

lucorum_count_exp02 <- sum(Guide$Experiment == "Experiment02" & Guide$Species == "lucorum", na.rm = TRUE)
print(lucorum_count_exp02)

cryptarum_count_exp02 <- sum(Guide$Experiment == "Experiment02" & Guide$Species == "cryptarum", na.rm = TRUE)
print(cryptarum_count_exp02)

terrestris_count_exp02 <- sum(Guide$Experiment == "Experiment02" & Guide$Species == "terrestris", na.rm = TRUE)
print(terrestris_count_exp02)

male_count_exp02 <- sum(Guide$Experiment == "Experiment02" & Guide$Male == 1, na.rm = TRUE)
print(male_count_exp02)

valid_mort_count_exp02 <- sum(Guide$Experiment == "Experiment02" & Guide$Valid_mort_data == 1, na.rm = TRUE)
print(valid_mort_count_exp02)

# Final treatment counts after all filtering.
print(table(Guide2$Treatment, useNA = "ifany"))

# These should align with the manuscript sample sizes for Experiment 2:
# ARCON = 20, WCON = 23, RAMI = 25, WAMI = 23,
# RCYN = 14, WCYN = 26, RDIM = 15, WDIM = 14.
# from the paper-facing mixed-effects model.

# Quick audit for missing core fields after filtering.
sum(is.na(Guide2$Time_in_experiment_hours))
sum(is.na(Guide2$Mortality))
sum(is.na(Guide2$Bee_weight_g))
sum(is.na(Guide2$Colony_ID))
sum(is.na(Guide2$Experimental_batch))

# Retain the legacy object name used in the original script.
GuideTer <- Guide2
head(GuideTer)
tail(GuideTer)

# Confirm that non-terrestris species have been excluded.
lucorum_count <- sum(GuideTer$Species == "lucorum", na.rm = TRUE)
print(lucorum_count)

cryptarum_count <- sum(GuideTer$Species == "cryptarum", na.rm = TRUE)
print(cryptarum_count)

magnus_count <- sum(GuideTer$Species == "magnus", na.rm = TRUE)
print(magnus_count)

pratorum_count <- sum(GuideTer$Species == "pratorum", na.rm = TRUE)
print(pratorum_count)

###############################################################################
# 2D. Survival inspection: Kaplan–Meier curves for the full Experiment 2 set
#
# These exploratory survival plots are used to confirm that treatment coding,
# filtering, and overall mortality patterns match the original script.
###############################################################################

# Quick survival plot for the full filtered Experiment 2 set

km_trt_fit_exp2 <- survfit(Surv(Time_in_experiment_hours, Mortality) ~ Treatment, data = GuideTer)
summary(km_trt_fit_exp2)

ggsurvplot(km_trt_fit_exp2, pval = TRUE, data = GuideTer)


###############################################################################
# 2E. Mixed-effects Cox modelling for Experiment 2
#
# Aim:
#   Model hazard of death over time while accounting for treatment, starting
#   bee weight, and non-independence among colonies / experimental batches.
#
# Important detail:
#   Wild bees do not have colony IDs in the same way as reared bees. To preserve
#   the original modelling logic and avoid dropping those rows, missing colony
#   IDs are coded as "Wild" before fitting the mixed-effects models.
###############################################################################

GuideTer$Colony_ID <- as.character(GuideTer$Colony_ID)
GuideTer$Colony_ID[is.na(GuideTer$Colony_ID)] <- "Wild"
GuideTer$Colony_ID <- factor(GuideTer$Colony_ID)
GuideTer$Experimental_batch <- factor(GuideTer$Experimental_batch)

# Full model including all filtered treatments.
FM <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + Bee_weight_g + (1 | as.factor(Colony_ID)) +
    (1 | as.factor(Experimental_batch)),
  data = GuideTer
)
summary(FM)
confint(FM)

# DIM is a positive control and was not retained in the manuscript-facing
# Experiment 2 model due to consistent 100% mortality
ddGuide <- droplevels(
  subset(
    GuideTer,
    Treatment != "WDIM" & Treatment != "RDIM"
  )
)
# dd = no DIM

###############################################################################
# 2F. Candidate model comparison
#
# Multiple biologically plausible models are fitted to evaluate whether bee
# weight and/or the random-effects structure materially improve the fit.
# The original workflow retained ddM2 as the simplest defensible model.
###############################################################################

ddFM <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + Bee_weight_g + (1 | as.factor(Colony_ID)) +
    (1 | as.factor(Experimental_batch)),
  data = ddGuide
)

ddM1 <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + (1 | as.factor(Colony_ID)) +
    (1 | as.factor(Experimental_batch)),
  data = ddGuide
)

ddM2 <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + Bee_weight_g +
    (1 | as.factor(Experimental_batch)),
  data = ddGuide
)

ddM3 <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + Bee_weight_g + (1 | as.factor(Colony_ID)),
  data = ddGuide
)

ddM4 <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + (1 | as.factor(Experimental_batch)),
  data = ddGuide
)

ddM5 <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + (1 | as.factor(Colony_ID)),
  data = ddGuide
)

ddM0 <- coxme(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    1 + Bee_weight_g + (1 | as.factor(Colony_ID)) +
    (1 | as.factor(Experimental_batch)),
  data = ddGuide
)

model.sel(ddFM, ddM1, ddM2, ddM3, ddM4, ddM5, ddM0)
# No single model is a dominant winner. The original workflow proceeds with
# ddM2 because it is the simplest model that still captures the main structure.

summary(ddM2)
confint(ddM2)

###############################################################################
# 2G. Planned contrasts for the manuscript table
#
# These are the six a priori comparisons used in the manuscript-facing results.
# Bonferroni correction is applied because the contrast set was predefined.
###############################################################################

ddM2emms <- emmeans(ddM2, ~ Treatment)

# All pairwise comparisons are available but are not the manuscript focus.
pairs(ddM2emms)

exp2_contrasts <- contrast(
  ddM2emms,
  method = list(
    "ARCON vs RAMI" = c(1, -1, 0, 0, 0, 0),
    "ARCON vs RCYN" = c(1, 0, -1, 0, 0, 0),
    "RAMI vs WAMI"  = c(0, 1, 0, -1, 0, 0),
    "RCYN vs WCYN"  = c(0, 0, 1, 0, -1, 0),
    "WCON vs WAMI"  = c(0, 0, 0, 1, 0, -1),
    "WCON vs WCYN"  = c(0, 0, 0, 0, 1, -1)
  ),
  adjust = "bonferroni"
)

summary(exp2_contrasts)
confint(exp2_contrasts)

###############################################################################
# 2H. Assumption checks for the retained mixed-effects model
#
# cox.zph() does not operate directly on coxme objects, so a fixed-effect coxph
# analogue is used for proportional hazards checks.
###############################################################################

ddM2_phcheck <- coxph(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~
    Treatment + Bee_weight_g + strata(Experimental_batch),
  data = ddGuide
)
cox.zph(ddM2_phcheck)
plot(cox.zph(ddM2_phcheck))
# The proportional hazards check is borderline, but this model remains the most
# defensible overall option for the manuscript-focused treatment contrasts.

# Linearity of the continuous covariate is checked using martingale residuals.
model_exp2_weight <- coxph(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~ Bee_weight_g,
  data = ddGuide
)
martingale_residuals <- residuals(model_exp2_weight, type = "martingale")
plot(
  ddGuide$Bee_weight_g,
  martingale_residuals,
  main = "Martingale Residuals vs. Bee_weight_g",
  xlab = "Bee_weight_g",
  ylab = "Martingale Residuals"
)
abline(h = 0, col = "red")

###############################################################################
# 2I. Figure 2
#
# The publication-facing figure retains DIM so the
# positive control remains visible. Significance letters are added manually to
# match the attached manuscript figure.
###############################################################################


noDelkm_trt_fit <- survfit(Surv(Time_in_experiment_hours, Mortality) ~ Treatment, data = GuideTer)
summary(noDelkm_trt_fit)

# Quick checks before building the publication figure.
ggsurvplot(noDelkm_trt_fit)

# Explicit plotting order for the manuscript-facing figure.
GuideTer$Treatment <- factor(
  GuideTer$Treatment,
  levels = c("ARCON", "RAMI", "RCYN", "RDIM", "WAMI", "WCON", "WCYN", "WDIM")
)

noDelkm_trt_fit <- survfit(
  Surv(Time_in_experiment_hours, Mortality) ~ Treatment,
  data = GuideTer
)

trt_cols_exp2 <- c(
  "ARCON" = "#0072B2",
  "RAMI"  = "#E69F00",
  "RCYN"  = "#56B4E9",
  "RDIM"  = "#000000",
  "WAMI"  = "#009E73",
  "WCON"  = "#F0E442",
  "WCYN"  = "#CC79A7",
  "WDIM"  = "#D55E00"
)

b <- ggsurvplot(
  noDelkm_trt_fit,
  data = GuideTer,
  size = 1.8,
  legend = "right",
  legend.labs = c(
    "Commercial-Control",
    "Commercial-Amistar",
    "Commercial-Cyantraniliprole",
    "Commercial-Positive Control",
    "Wild-Amistar",
    "Wild-Control",
    "Wild-Cyantraniliprole",
    "Wild-Positive Control"
  ),
  legend.title = "",
  font.legend = 18,
  font.x = 26,
  font.y = 26,
  font.xtickslab = 18,
  font.ytickslab = 18,
  break.time.by = 24,
  xlim = c(0, 50),
  ylim = c(0, 1.05),
  censor = FALSE,
  conf.int = FALSE,
  xlab = "Time (hours)",
  ylab = "Survival",
  surv.scale = "percent",
  palette = unname(trt_cols_exp2[c("ARCON", "RAMI", "RCYN", "RDIM", "WAMI", "WCON", "WCYN", "WDIM")]),
  ggtheme = theme_classic(base_size = 18) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(colour = "black"),
      legend.position = "right",
      legend.text = element_text(size = 18, colour = "black"),
      legend.key = element_blank(),
      plot.margin = margin(t = 10, r = 65, b = 10, l = 10)
    )
)

# Significance letters at 48 h, as shown in the manuscript figure.
km_48_exp2 <- summary(noDelkm_trt_fit, times = 48, extend = TRUE)

letters_df_exp2 <- data.frame(
  Treatment = gsub("^Treatment=", "", km_48_exp2$strata),
  survival = km_48_exp2$surv,
  stringsAsFactors = FALSE
)

letters_df_exp2 <- letters_df_exp2[
  match(c("ARCON", "RAMI", "RCYN", "RDIM", "WAMI", "WCON", "WCYN", "WDIM"),
        letters_df_exp2$Treatment),
]

letters_df_exp2$label <- c("a", "b", "b", "c", "a", "a", "a", "a")
letters_df_exp2$x <- c(49.3, 49.3, 49.3, 5.2, 49.3, 49.3, 49.3, 49.3)
letters_df_exp2$y <- c(
  letters_df_exp2$survival[letters_df_exp2$Treatment == "ARCON"] + 0.005,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "RAMI"]  + 0.007,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "RCYN"]  + 0.004,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "RDIM"]  + 0.015,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "WAMI"]  - 0.002,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "WCON"]  - 0.001,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "WCYN"]  - 0.002,
  letters_df_exp2$survival[letters_df_exp2$Treatment == "WDIM"]  - 0.002
)
letters_df_exp2$colour <- trt_cols_exp2[letters_df_exp2$Treatment]

for (i in seq_len(nrow(letters_df_exp2))) {
  b$plot <- b$plot +
    annotate(
      "text",
      x = letters_df_exp2$x[i],
      y = letters_df_exp2$y[i],
      label = letters_df_exp2$label[i],
      colour = letters_df_exp2$colour[i],
      size = 8,
      hjust = 0,
      family = "sans"
    )
}

print(b$plot)

ggsave(
  filename = "wild_vs_reared_bt_worker_mort_exp2.pdf",
  plot = b$plot,
  width = 11,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "wild_vs_reared_bt_worker_mort_exp2.png",
  plot = b$plot,
  width = 11,
  height = 7.5,
  dpi = 300,
  bg = "white"
)


###############################################################################
# EXPERIMENT 3
# Outdoor-reared vs lab-reared commercial worker bumble bees
# Goal: test whether a period of outdoor exposure ('wild-hardening') reduces
# the elevated pesticide susceptibility typically seen in lab-reared commercial
# bees.
#
# General logic for this experiment:
#   1. Filter the consolidated dataset to the Experiment 3 analysis subset
#   2. Audit sample sizes and key variables after filtering
#   3. Inspect survival curves over time
#   4. Attempt a Cox model, but retain chi-squared / Fisher testing as the main
#      manuscript-facing mortality analysis because mortality patterns are too
#      sparse / separated for a stable proportional-hazards approach
#   5. Recreate manuscript Figure 3 with manual significance lettering
#   6. Explore whether bee weight helps explain the mortality pattern
###############################################################################

###############################################################################
# 3A. Plotting palette for Experiment 3
#
# These colours are retained from the working script so the legend order and
# survival curves remain visually consistent with the manuscript figure.
###############################################################################

cbPaletteOkabeIto_exp3 <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000")

###############################################################################
# 3B. Check that the consolidated dataset contains the columns needed for
#     Experiment 3
###############################################################################

required_cols_exp3 <- c(
  "Experiment", "Experimental_batch", "Caste", "Treatment", "Species", "Male",
  "Bee_weight_g", "Time_in_experiment_hours", "Mortality", "Valid_mort_data"
)
missing_cols_exp3 <- setdiff(required_cols_exp3, names(Guide))
if (length(missing_cols_exp3) > 0) {
  stop(paste("The following required columns are missing for Experiment 3:",
             paste(missing_cols_exp3, collapse = ", ")))
}

###############################################################################
# 3C. Filter the consolidated dataset to the Experiment 3 analysis set
#
# Paper-wide inclusion rules applied here:
#   - Experiment03 only
#   - Worker caste only
#   - Bombus terrestris only
#   - Males excluded
#   - Valid mortality data only
###############################################################################

Guide3 <- subset(
  Guide,
  Experiment == "Experiment03" &
    Caste == "Worker" &
    Species == "terrestris" &
    Male != 1 &
    Valid_mort_data == 1
)

# Force treatment order so tables, legends, and manual significance labels all
# refer to treatments in a stable, known order.
Guide3$Treatment <- factor(
  Guide3$Treatment,
  levels = c("LBC", "LBD", "LBT", "ODC", "ODD", "ODT")
)

###############################################################################
# 3D. Audit the Experiment 3 analysis subset

###############################################################################

exp03_count <- sum(Guide$Experiment == "Experiment03", na.rm = TRUE)
print(exp03_count)

worker_count_exp03 <- sum(Guide$Experiment == "Experiment03" & Guide$Caste == "Worker", na.rm = TRUE)
print(worker_count_exp03)

lucorum_count_exp03 <- sum(Guide$Experiment == "Experiment03" & Guide$Species == "lucorum", na.rm = TRUE)
print(lucorum_count_exp03)

terrestris_count_exp03 <- sum(Guide$Experiment == "Experiment03" & Guide$Species == "terrestris", na.rm = TRUE)
print(terrestris_count_exp03)

male_count_exp03 <- sum(Guide$Experiment == "Experiment03" & Guide$Male == 1, na.rm = TRUE)
print(male_count_exp03)

valid_mort_count_exp03 <- sum(Guide$Experiment == "Experiment03" & Guide$Valid_mort_data == 1, na.rm = TRUE)
print(valid_mort_count_exp03)

# Final treatment counts after all filtering
print(table(Guide3$Treatment, useNA = "ifany"))

# These should match the manuscript/sample sizes for Experiment 3:
# LBC = 43, LBT = 38, LBD = 23, ODC = 32, ODT = 32, ODD = 17

# Quick audit for missing core model / plotting fields after filtering
sum(is.na(Guide3$Time_in_experiment_hours))
sum(is.na(Guide3$Mortality))
sum(is.na(Guide3$Bee_weight_g))

###############################################################################
# 3E. Survival inspection
#
# These Kaplan–Meier curves are used first as a visual check that the filtered
# data show the expected mortality pattern. They are exploratory checks rather
# than the final manuscript figure.
###############################################################################

km_trt_fit_exp3 <- survfit(Surv(Time_in_experiment_hours, Mortality) ~ Treatment, data = Guide3)
summary(km_trt_fit_exp3)

# Quick visual checks retained from the working script.
ggsurvplot(km_trt_fit_exp3, data = Guide3)
ggsurvplot(km_trt_fit_exp3, pval = TRUE, data = Guide3)

###############################################################################
# 3F. Mortality inference for the manuscript
#
# Positive controls are removed from the inferential comparison set because they
# are included to verify the assay response, not to address the main biological
# treatment questions.
###############################################################################

NoDimGuide3 <- droplevels(subset(Guide3, Treatment != "ODD" & Treatment != "LBD"))

# A Cox model is attempted first because it preserves time-to-event information.
# This was judged unsuitable / unstable for the final
# manuscript analysis because the mortality pattern is too sparse and separated.
PH3 <- coxph(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~ as.factor(Treatment) +
    Bee_weight_g,
  data = NoDimGuide3
)
summary(PH3)

# Main manuscript-facing approach:
# chi-squared tables (with Fisher tests where counts are sparse) because these
# data are strongly bifurcated into high-mortality vs near-zero-mortality groups.
NoDIMtbl3 <- table(NoDimGuide3$Treatment, NoDimGuide3$Mortality)
NoDIMtbl3
# 0 = survived
# 1 = died

NoDIMcs3 <- chisq.test(NoDIMtbl3)
NoDIMcs3
# Global treatment effect is expected to be highly significant.

# There are six pairwise comparisons in Experiment 3.
# Bonferroni-adjusted alpha: 0.05 / 6 = 0.008333.

# Helper function:
# After subsetting a factor, unused levels can remain attached. droplevels()
# removes those unused levels so the resulting contingency table contains only
# the two groups being compared.
make_pair_table_exp3 <- function(data, grp1, grp2) {
  subdat <- droplevels(subset(data, Treatment %in% c(grp1, grp2)))
  table(subdat$Treatment, subdat$Mortality)
}

# 1. LBC vs LBT
LBCLBTtbl <- make_pair_table_exp3(Guide3, "LBC", "LBT")
LBCLBTtbl
LBCLBTcs <- chisq.test(LBCLBTtbl)
LBCLBTcs
# This should be significant: lab control vs lab pesticide.

# 2. LBC vs ODC
LBCODCtbl <- make_pair_table_exp3(Guide3, "LBC", "ODC")
LBCODCtbl
LBCODCcs <- chisq.test(LBCODCtbl)
LBCODCcs
# This comparison is sparse enough that Fisher's exact test is also useful.
LBCODCfisher_result <- fisher.test(LBCODCtbl)
print(LBCODCfisher_result)
# In practice this comparison is treated as non-significant.

# 3. LBC vs ODT
LBCODTtbl <- make_pair_table_exp3(Guide3, "LBC", "ODT")
LBCODTtbl
LBCODTcs <- chisq.test(LBCODTtbl)
LBCODTcs

# 4. ODC vs ODT
ODCODTtbl <- make_pair_table_exp3(Guide3, "ODC", "ODT")
ODCODTtbl
ODCODTcs <- chisq.test(ODCODTtbl)
ODCODTcs

# 5. ODC vs LBT
ODCLBTtbl <- make_pair_table_exp3(Guide3, "ODC", "LBT")
ODCLBTtbl
ODCLBTcs <- chisq.test(ODCLBTtbl)
ODCLBTcs
# This comparison stretches assumptions somewhat, but the effect is strong.

# 6. LBT vs ODT
LBTODTtbl <- make_pair_table_exp3(Guide3, "LBT", "ODT")
LBTODTtbl
LBTODTcs <- chisq.test(LBTODTtbl)
LBTODTcs
# This comparison is non-significant by chi-squared testing.

LBTODTfisher_result <- fisher.test(LBTODTtbl)
print(LBTODTfisher_result)
# Fisher's test is significant, but we retain the chi-squared
# interpretation for consistency with the main workflow, as the 
# chi-square model di fit

# Summary of pairwise outputs retained for quick review against manuscript notes
LBCLBTcs
LBCODCcs
print(LBCODCfisher_result)
LBCODTcs
ODCODTcs
ODCLBTcs
LBTODTcs
print(LBTODTfisher_result)

###############################################################################
# 3G. Figure 3
#
# Purpose:
# Recreate the manuscript survival figure for Experiment 3.
#
# Notes:
# - The standard survival plot is retained from the working script.
# - Significance letters are added manually at 48 h to match the figure used in
#   the manuscript.
# - Letter positions are hand-tuned and may need small adjustment if export size
#   is changed.
###############################################################################

a3 <- ggsurvplot(
  data = Guide3,
  km_trt_fit_exp3,
  size = 1.5,
  legend = "right",
  legend.labs = c(
    "Lab - Commercial - Control",
    "Lab - Commercial - Positive Control",
    "Lab - Commercial - Pesticide",
    "Outdoor - Commercial - Control",
    "Outdoor - Commercial - Positive Control",
    "Outdoor - Commercial - Pesticide"
  ),
  font.legend = 15,
  legend.title = c(""),
  font.x = c(20),
  font.y = c(20),
  font.xtickslab = c(15),
  break.x.by = 24,
  font.ytickslab = c(15),
  censor.shape = "",
  xlab = ("Time (hours)"),
  ylab = ("Survival (%)"),
  palette = cbPaletteOkabeIto_exp3,
  ggtheme = theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(hjust = 0.4)
    ) +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black")
    )
)
a3
# Colours were adjusted in post-production in the original workflow, but the
# plotted survival content is the same.

# Significance letters, as shown in the manuscript figure:
# Lab-Control = a
# Lab-Positive Control = b
# Lab-Pesticide = b
# Outdoor-Control = a
# Outdoor-Positive Control = b
# Outdoor-Pesticide = a
km_48_exp3 <- summary(km_trt_fit_exp3, times = 48, extend = TRUE)

letters_df_exp3 <- data.frame(
  Treatment = gsub("^Treatment=", "", km_48_exp3$strata),
  survival = km_48_exp3$surv,
  stringsAsFactors = FALSE
)

letters_df_exp3 <- letters_df_exp3[match(c("LBC", "LBD", "LBT", "ODC", "ODD", "ODT"),
                                         letters_df_exp3$Treatment), ]
letters_df_exp3$label <- c("a", "b", "b", "a", "b", "a")
letters_df_exp3$x <- c(49.3, 5.2, 49.3, 49.3, 5.2, 49.3)
letters_df_exp3$y <- c(
  letters_df_exp3$survival[letters_df_exp3$Treatment == "LBC"] + 0.005,
  letters_df_exp3$survival[letters_df_exp3$Treatment == "LBD"] + 0.015,
  letters_df_exp3$survival[letters_df_exp3$Treatment == "LBT"] + 0.007,
  letters_df_exp3$survival[letters_df_exp3$Treatment == "ODC"] - 0.001,
  letters_df_exp3$survival[letters_df_exp3$Treatment == "ODD"] + 0.015,
  letters_df_exp3$survival[letters_df_exp3$Treatment == "ODT"] - 0.002
)

letters_df_exp3$colour <- c(
  "LBC" = "#0072B2",
  "LBD" = "#E69F00",
  "LBT" = "#56B4E9",
  "ODC" = "#009E73",
  "ODD" = "#F0E442",
  "ODT" = "#000000"
)[letters_df_exp3$Treatment]

for (i in seq_len(nrow(letters_df_exp3))) {
  a3$plot <- a3$plot +
    annotate(
      "text",
      x = letters_df_exp3$x[i],
      y = letters_df_exp3$y[i],
      label = letters_df_exp3$label[i],
      colour = letters_df_exp3$colour[i],
      size = 8,
      hjust = 0,
      family = "sans"
    )
}

print(a3$plot)

ggsave("Experiment03_lab_into_wild_worker_mortality.pdf",
       plot = a3$plot, width = 11, height = 7.5, dpi = 300, bg = "white")
ggsave("Experiment03_lab_into_wild_worker_mortality.png",
       plot = a3$plot, width = 11, height = 7.5, dpi = 300, bg = "white")

###############################################################################
# 3H. Additional explanatory-factor checks: weight
#
# These analyses are retained from the original script because they address the
# alternative explanation that mortality differences are primarily driven by
# body size rather than treatment / rearing history.
###############################################################################

# Overall relationship between treatment, mortality, and body weight
plot(Bee_weight_g ~ Mortality, data = NoDimGuide3)
plot(Bee_weight_g ~ as.factor(Treatment), data = NoDimGuide3)
# Outdoor-reared bees appear lighter overall.

WeightComp3 <- aov(Bee_weight_g ~ as.factor(Treatment),
                   data = NoDimGuide3)
summary(WeightComp3)
WeightCompTukey_Results3 <- TukeyHSD(WeightComp3)
print(WeightCompTukey_Results3)
# Outdoor-reared groups are lighter than lab-reared groups.

# Mortality across the full experiment is too concentrated in one treatment to
# model covariate effects sensibly, so the weight check is restricted to the two
# pesticide-exposed groups.
PesticideOnlyGuide3 <- droplevels(subset(NoDimGuide3, Treatment == "LBT" | Treatment == "ODT"))

plot(Bee_weight_g ~ Mortality, data = PesticideOnlyGuide3)
plot(Bee_weight_g ~ as.factor(Treatment), data = PesticideOnlyGuide3)

PesticideMortVSWeight3 <- glm(Bee_weight_g ~ as.factor(Mortality),
                              data = PesticideOnlyGuide3)
summary(PesticideMortVSWeight3)
confint(PesticideMortVSWeight3)
# Non-significant.

# Scatter plot of body weight vs survival status for the two pesticide-exposed
# groups, retained as a visual check.
LBTODTGuide <- droplevels(subset(Guide3, Treatment == "LBT" | Treatment == "ODT"))
head(LBTODTGuide)

f3 <- ggplot(data = LBTODTGuide, aes(y = as.factor(Mortality), x = Bee_weight_g, color = Treatment)) +
  geom_point(position = position_jitter(width = 0, height = 0.1), size = 2) +
  theme_bw() +
  coord_flip() +
  scale_x_continuous(limits = c(0, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = "right",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)
  ) +
  labs(x = "Bee Weight at Start (g)", y = "Survival") +
  scale_color_manual(values = c("ODT" = "#56B4E9", "LBT" = "#000000"))
f3

g3 <- f3 + scale_y_discrete(labels = c("0" = "Survived", "1" = "Died"))
g3
ggsave("Experiment03_lab_into_wild_worker_weight_LBTODT_scatter.pdf", g3, width = 6, height = 8, dpi = 300)


###############################################################################
# EXPERIMENT 4
# Lab-reared wild-genotype vs. lab-reared commercial worker bumble bees
#
# Aim:
# Test whether wild-genotype bees become more susceptible when reared under lab
# conditions, using the same general framework as Experiment 3.
#
# Manuscript logic:
# - Survival curves are used to inspect mortality trajectories over time.
# - Positive controls are retained in figures for context.
# - Positive controls are removed from manuscript-facing mortality comparisons.
# - A Cox model is explored, but chi-squared / Fisher tests remain the main
#   inferential workflow because the original analysis was built around overall
#   mortality outcomes rather than detailed time-to-event modelling.
###############################################################################

###############################################################################
# 4A. Filter consolidated dataset to the Experiment 4 analysis set
###############################################################################

# Colour palette used for the Experiment 4 figure. Colours follow the same logic
# as the original script and are retained for consistency with the attached
# publication figure.
cbPaletteOkabeIto_exp4 <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000")

# Check that all columns needed for the Experiment 4 workflow are present.
required_cols_exp4 <- c(
  "Experiment", "Experimental_batch", "Caste", "Treatment", "Species", "Male",
  "Bee_weight_g", "Time_in_experiment_hours", "Mortality", "Valid_mort_data"
)
missing_cols_exp4 <- setdiff(required_cols_exp4, names(Guide))
if (length(missing_cols_exp4) > 0) {
  stop(paste("The following required columns are missing for Experiment 4:",
             paste(missing_cols_exp4, collapse = ", ")))
}

# Construct the Experiment 4 analysis dataset from the consolidated CSV using
# the paper-wide inclusion rules.
Guide4 <- subset(
  Guide,
  Experiment == "Experiment04" &
    Caste == "Worker" &
    Species == "terrestris" &
    Male != 1 &
    Valid_mort_data == 1
)

# Force treatment order so legends, tables, and figure annotations remain
# stable and match the intended experiment logic.
Guide4$Treatment <- factor(
  Guide4$Treatment,
  levels = c("KC", "KD", "KT", "WC", "WD", "WT")
)

###############################################################################
# 4B. Audit sample sizes and key variables after filtering
###############################################################################

# These checks are retained from the original workflow as a sanity check that
# the filtered Experiment 4 dataset matches the manuscript analysis set.
exp04_count <- sum(Guide$Experiment == "Experiment04", na.rm = TRUE)
print(exp04_count)

worker_count_exp04 <- sum(Guide$Experiment == "Experiment04" & Guide$Caste == "Worker", na.rm = TRUE)
print(worker_count_exp04)

lucorum_count_exp04 <- sum(Guide$Experiment == "Experiment04" & Guide$Species == "lucorum", na.rm = TRUE)
print(lucorum_count_exp04)

terrestris_count_exp04 <- sum(Guide$Experiment == "Experiment04" & Guide$Species == "terrestris", na.rm = TRUE)
print(terrestris_count_exp04)

male_count_exp04 <- sum(Guide$Experiment == "Experiment04" & Guide$Male == 1, na.rm = TRUE)
print(male_count_exp04)

valid_mort_count_exp04 <- sum(Guide$Experiment == "Experiment04" & Guide$Valid_mort_data == 1, na.rm = TRUE)
print(valid_mort_count_exp04)

# Final treatment counts after all filtering.
print(table(Guide4$Treatment, useNA = "ifany"))

# These should match the manuscript sample sizes for Experiment 4:
# KC = 27, KT = 30, KD = 23, WC = 35, WT = 38, WD = 27

# Quick audit for missing values in the key mortality-analysis fields.
sum(is.na(Guide4$Time_in_experiment_hours))
sum(is.na(Guide4$Mortality))
sum(is.na(Guide4$Bee_weight_g))

###############################################################################
# 4C. Survival inspection: Kaplan–Meier curves
#
# These plots are retained as the first visual check that the expected survival
# pattern is present before the manuscript-facing inferential tests are run.
###############################################################################

km_trt_fit_exp4 <- survfit(Surv(Time_in_experiment_hours, Mortality) ~ Treatment, data = Guide4)
summary(km_trt_fit_exp4)

# Quick visual checks mirroring the original script. These are exploratory and
# are not manuscript figures.
ggsurvplot(km_trt_fit_exp4, data = Guide4)

###############################################################################
# 4D. Mortality inference for manuscript table
###############################################################################

# Positive controls are excluded from the manuscript-facing mortality
# comparisons because they are included to verify assay sensitivity rather than
# to address the focal biological comparisons.
NoDimGuide4 <- droplevels(subset(Guide4, Treatment != "KD" & Treatment != "WD"))

# As in Experiment 3, a Cox model is explored first as a diagnostic check.
# This preserves the logic of the original script and confirms whether a more
# time-explicit approach is plausible.
PH4 <- coxph(
  Surv(time = Time_in_experiment_hours, event = Mortality) ~ as.factor(Treatment) +
    Bee_weight_g,
  data = NoDimGuide4
)
summary(PH4)
# Unlike Experiment 3, this model converges. The chi-squared workflow is still
# retained because that is the original manuscript-facing analysis and improves
# comparability across experiments.

# Global treatment effect on overall mortality, excluding positive controls.
NoDIMtbl4 <- table(NoDimGuide4$Treatment, NoDimGuide4$Mortality)
NoDIMtbl4
# 0 = survived / censored, 1 = died

NoDIMcs4 <- chisq.test(NoDIMtbl4)
NoDIMcs4
# Strong global treatment effect.

# Six pairwise treatment comparisons were planned. The Bonferroni-adjusted
# threshold is therefore 0.05 / 6 = 0.008333.
#
#   1   WC vs WT
#   2   WC vs KT
#   3   WC vs KC
#   4   KC vs KT
#   5   KC vs WT
#   6   KT vs WT

make_pair_table_exp4 <- function(data, grp1, grp2) {
  subdat <- droplevels(subset(data, Treatment %in% c(grp1, grp2)))
  table(subdat$Treatment, subdat$Mortality)
}

# 1. WC vs WT
WCWTtbl <- make_pair_table_exp4(Guide4, "WC", "WT")
WCWTtbl
WCWTcs <- chisq.test(WCWTtbl)
WCWTcs
# Expected significant control vs pesticide difference.

# 2. WC vs KT
WCKTtbl <- make_pair_table_exp4(Guide4, "WC", "KT")
WCKTtbl
WCKTcs <- chisq.test(WCKTtbl)
WCKTcs
# Expected significant control vs pesticide difference.

# 3. WC vs KC
WCKCtbl <- make_pair_table_exp4(Guide4, "WC", "KC")
WCKCtbl
WCKCcs <- chisq.test(WCKCtbl)
WCKCcs
# This comparison is effectively control vs control and is expected to be non-significant.
# Fisher's exact test is retained because the sparse structure can weaken the chi-squared assumptions.
WCKCfisher_result <- fisher.test(WCKCtbl)
print(WCKCfisher_result)

# 4. KC vs KT
KCKTtbl <- make_pair_table_exp4(Guide4, "KC", "KT")
KCKTtbl
KCKTcs <- chisq.test(KCKTtbl)
KCKTcs
# Expected significant commercial control vs pesticide difference.

# 5. KC vs WT
KCWTtbl <- make_pair_table_exp4(Guide4, "KC", "WT")
KCWTtbl
KCWTcs <- chisq.test(KCWTtbl)
KCWTcs
# Expected significant control vs pesticide difference.

# 6. KT vs WT
KTWTtbl <- make_pair_table_exp4(Guide4, "KT", "WT")
KTWTtbl
KTWTcs <- chisq.test(KTWTtbl)
KTWTcs
# Borderline at the raw p-value level, but not significant after the planned
# Bonferroni-adjusted threshold.

# Compact summary block for quick inspection.
WCWTcs
WCKTcs
WCKCcs
print(WCKCfisher_result)
KCKTcs
KCWTcs
KTWTcs

###############################################################################
# 4E. Publication figure (Figure 4)
#
# Purpose:
# Create the manuscript survival figure for Experiment 4.
#
# Notes:
# - Positive controls are retained in the survival plot for context.
# - Significance letters are added manually to match the attached figure.
# - Only the control and focal pesticide groups carry visible letters in the
#   figure body; the positive controls are represented in the legend.
###############################################################################

a4 <- ggsurvplot(
  data = Guide4,
  km_trt_fit_exp4,
  size = 1.5,
  legend = "right",
  legend.labs = c(
    "Lab - Commercial - Control",
    "Lab - Commercial - Positive Control",
    "Lab - Commercial - Pesticide",
    "Lab - Wild - Control",
    "Lab - Wild - Positive Control",
    "Lab - Wild - Pesticide"
  ),
  font.legend = 15,
  legend.title = c(""),
  font.x = c(20),
  font.y = c(20),
  font.xtickslab = c(15),
  break.x.by = 24,
  font.ytickslab = c(15),
  censor.shape = "",
  xlab = "Time (hours)",
  ylab = "Survival (%)",
  palette = cbPaletteOkabeIto_exp4,
  ggtheme = theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.title.x = element_text(hjust = 0.4)
    ) +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black")
    )
)

# Manual significance letters to match the attached figure.
# Commercial-Control = a
# Commercial-Pesticide = b
# Wild-Control = a
# Wild-Pesticide = b
km_48_exp4 <- summary(km_trt_fit_exp4, times = 48, extend = TRUE)

letters_df_exp4 <- data.frame(
  Treatment = gsub("^Treatment=", "", km_48_exp4$strata),
  survival = km_48_exp4$surv,
  stringsAsFactors = FALSE
)

letters_df_exp4 <- letters_df_exp4[match(c("KC", "KD", "KT", "WC", "WD", "WT"),
                                         letters_df_exp4$Treatment), ]

# Only the focal control and pesticide treatments receive visible letters.
letters_df_exp4$label <- c("a", NA, "b", "a", NA, "b")

# Text is placed on the native 0-1 survival scale because ggsurvplot handles
# the percent display internally.
letters_df_exp4$x <- c(49.3, 5.2, 49.3, 49.3, 5.2, 49.3)
letters_df_exp4$y <- c(
  letters_df_exp4$survival[letters_df_exp4$Treatment == "KC"] + 0.005,
  letters_df_exp4$survival[letters_df_exp4$Treatment == "KD"] + 0.015,
  letters_df_exp4$survival[letters_df_exp4$Treatment == "KT"] + 0.007,
  letters_df_exp4$survival[letters_df_exp4$Treatment == "WC"] - 0.001,
  letters_df_exp4$survival[letters_df_exp4$Treatment == "WD"] + 0.015,
  letters_df_exp4$survival[letters_df_exp4$Treatment == "WT"] - 0.002
)

letters_df_exp4$colour <- c(
  "KC" = "#0072B2",
  "KD" = "#E69F00",
  "KT" = "#56B4E9",
  "WC" = "#009E73",
  "WD" = "#F0E442",
  "WT" = "#000000"
)[letters_df_exp4$Treatment]

for (i in seq_len(nrow(letters_df_exp4))) {
  if (!is.na(letters_df_exp4$label[i])) {
    a4$plot <- a4$plot +
      annotate(
        "text",
        x = letters_df_exp4$x[i],
        y = letters_df_exp4$y[i],
        label = letters_df_exp4$label[i],
        colour = letters_df_exp4$colour[i],
        size = 8,
        hjust = 0,
        family = "sans"
      )
  }
}

print(a4$plot)

ggsave("Experiment04_wild_into_lab_worker_mortality.pdf",
       plot = a4$plot, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("Experiment04_wild_into_lab_worker_mortality.png",
       plot = a4$plot, width = 12, height = 8, dpi = 300, bg = "white")

###############################################################################
# 4F. Additional explanatory-factor checks: weight
#
# These retained checks ask whether Experiment 4 mortality patterns might be
# explained by body size differences rather than treatment / rearing history.
###############################################################################

# Overall relationship between treatment, mortality, and body weight.
plot(Bee_weight_g ~ Mortality, data = NoDimGuide4)
plot(Bee_weight_g ~ as.factor(Treatment), data = NoDimGuide4)
# Wild-genotype groups appear heavier overall.

WeightComp4 <- aov(Bee_weight_g ~ as.factor(Treatment),
                   data = NoDimGuide4)
summary(WeightComp4)
WeightCompTukey_Results4 <- TukeyHSD(WeightComp4)
print(WeightCompTukey_Results4)
# Significant differences in mean body weight among treatment groups.

# Because mortality is concentrated in the pesticide groups, the weight-vs-
# mortality check is restricted to KT and WT.
PesticideOnlyGuide4 <- droplevels(subset(NoDimGuide4, Treatment == "KT" | Treatment == "WT"))

plot(Bee_weight_g ~ Mortality, data = PesticideOnlyGuide4)
plot(Bee_weight_g ~ as.factor(Treatment), data = PesticideOnlyGuide4)

PesticideMortVSWeight4 <- glm(Bee_weight_g ~ as.factor(Mortality),
                              data = PesticideOnlyGuide4)
summary(PesticideMortVSWeight4)
confint(PesticideMortVSWeight4)
# This effect was non-significant.

# Scatter plot of body weight vs survival status for the two pesticide-exposed
# groups, retained as a visual check.
KTWTGuide4 <- droplevels(subset(Guide4, Treatment == "KT" | Treatment == "WT"))
head(KTWTGuide4)

f4 <- ggplot(data = KTWTGuide4, aes(y = as.factor(Mortality), x = Bee_weight_g, color = Treatment)) +
  geom_point(position = position_jitter(width = 0, height = 0.1), size = 2) +
  theme_bw() +
  coord_flip() +
  scale_x_continuous(limits = c(0, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = "right",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)
  ) +
  labs(x = "Bee Weight at Start (g)", y = "Survival") +
  scale_color_manual(values = c("WT" = "#56B4E9", "KT" = "#000000"))
f4

g4 <- f4 + scale_y_discrete(labels = c("0" = "Survived", "1" = "Died"))
g4
ggsave("Experiment04_wild_into_lab_worker_weight_KTWT_scatter.pdf", g4, width = 6, height = 8, dpi = 300)

###############################################################################
