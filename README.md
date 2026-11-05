# Hybrid-VR vs real-world object lifting
 
Code for:
 
> Sager CA, Zenti J, Marneweck M. A behavioral comparison of hybrid virtual reality and real-world object lifting. *Journal of Neurophysiology* (under review).
 
Motor Skill Lab, University of Oregon.
 
## Contents
 
| Script | What it does | Manuscript |
|---|---|---|
| `src/collection/rep_vr.py` | Runs the asymmetric-load lifting task (Vizard); records forces and object motion | Methods: Experimental Procedures |
| `src/extraction/behaviorR01.py` | Extracts lift-onset grip/lift forces, COP and M_Com from raw force, motion and trigger files; writes one summary file per participant | Methods: Data Processing |
| `src/extraction/vr_single_analysis.m` | Separates pre-switch and post-switch trials for one participant | Methods: Data Processing |
| `src/extraction/vr_group_analysis.m` | Runs `vr_single_analysis` for all participants | Methods: Data Processing |
| `src/stats/repVR_nonlinear_model_v3.R` | Nonlinear mixed-effects model of pre-switch M_Com | Results 3.1; Supplemental Tables S1–S2 |
| `src/figs/pre_5_row_single.m` | Separates pre-switch and post-switch trials for one participant | Methods: Data Processing, Figure 2 |
| `src/figs/pre_5_row_group.m` | Separates pre-switch and post-switch trials for all participant | Methods: Data Processing, Figure 2 |
| `src/stats/components.R` | Repeated-measures ANOVAs for COP_diff, LF_diff, GF_mean and R_COP | Results 3.1; Figure 2; Supplemental Figure S1 |
| `src/stats/correlation_co_r.R` | Group position-force correlations (Fisher z), 95% CIs and environment comparison | Results 3.2 |
| `src/stats/post_switch_ANOVA.R` | Mixed ANOVA on post-switch M_Com | Results 3.3; Figure 4 |
| `src/figs/individual_scatterplots_combined_v2.m` | Individual-participant COP_diff vs LF_diff scatterplots and per-participant r | Figure 3A |
| `src/figs/group_scatterplot_cop_lf_v2.m` | Group COP_diff vs LF_diff scatterplots with fitted line and 95% CI | Figure 3B–C |
 
## Software
 
- R 4.6.0 with nlme 3.1-169, afex 1.5-1, emmeans 2.0.3, openxlsx and tidyverse
- Python 3.12.12 with NumPy and Matplotlib; Vizard (data collection)
- MATLAB 26.1.0.3276743 (R2026a) Update 3
- Unity 2018.4.10f1 with SteamVR (hybrid-VR environment; not included)
 
## Running the code
 
**Data collection.** `rep_vr.py` requires Vizard, the PPT-N motion tracking system, and the force/torque sensors. At start-up, enter the subject number (with leading zero), trial start number, and run number; trial order is read from the condition files.
 
**Extraction.** Run `behaviorR01.py` for each participant to create the summary files, then run `vr_group_analysis.m` from the folder containing the summary files to extract pre-switch and post-switch M_Com for all participants.
 
**Statistics.** The R scripts read Excel files with one column per participant and one row per trial (pre-switch: 25 rows per environment, Trials 1–5 within Blocks 1–5) or per block (post-switch). Set the folder path at the top of each script. Instructions for running the nonlinear mixed-effects model are at the top of `repVR_nonlinear_model_v3.R`.
 
**Figures.** Run the scatterplot scripts from the folder containing the participant summary files. Within each script, set the marker to `'filled'` for real-world or `'o'` for hybrid-VR, and select the trial subset (all trials, last three trials of each block, or last three blocks).
 
## Data availability
 
Data are available from the corresponding author (Michelle Marneweck, mmar@uoregon.edu) upon reasonable request.
 
## License
 
Apache License 2.0 (see `LICENSE`).
