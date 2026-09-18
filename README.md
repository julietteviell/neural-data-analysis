# Classification

This MATLAB script processes electrophysiological recordings from dorsal horn (DH) spinal cord neurons, classifies cells based on their responsiveness to sensory stimuli, and plots peri-event activity (raw/normalized and z-scored) for ipsilateral vs. contralateral paw stimulation, split by cell category and by event type.

## Overview of the workflow

1. **Data loading & organization** — loads per-mouse recording files, splits units into left (`L`) and right (`R`) dorsal horn groups based on channel numbering, and reorganizes them into ipsi/contra (`LDH`/`RDH`) structures.
2. **Binned spike counts** — computes peri-event time histograms (raw and normalized) around each stimulus/withdrawal event, for both paws.
3. **Concatenation across animals** — merges per-animal binned data into single population-level structures (`miceall`, `normiceall`).
4. **Cell classification** — z-scores each neuron's response, classifies neurons as activated/inhibited/unresponsive per event, and builds ordered response heatmaps.
5. **Plotting** — generates two summary figures (z-scored and normalized activity) with **cell categories as columns** and **event types as rows**, comparing ipsilateral vs. contralateral paw responses.

## Functions used

### Local / defined in this script

| Function | Purpose |
|---|---|
| `computeMeanSEM(y)` | Computes the mean and standard error of the mean (SEM) across rows (neurons/trials) of a data matrix `y`. Returns the raw row if only one row is present (single-unit case), otherwise averages across the first dimension and divides the standard deviation by `sqrt(N)`. Used to build the mean traces and error bands plotted with `shadedErrorBar`. |

### External / dependency functions (defined elsewhere in the analysis pipeline)

| Function | Purpose |
|---|---|
| `crosco(file, p)` | Removes duplicate units across recordings based on co-firing probability, returning a cleaned unit structure `x`. |
| `freqVarRAMalg(var3d, varstring, b, bstring, p, file, edges, bin)` | Computes binned spike counts (mean firing rate) around each behavioral/stimulus event, returning both raw (`total_counts_mean_var`) and normalized (`total_counts_mean_var_norm`) peri-event time histograms per neuron. |
| `zscoreVarSCRAMalgo(varstring, bstring, mice, g)` | Z-scores each neuron's peri-event response relative to its own baseline activity; also returns the baseline (`basal`) values used for normalization. |
| `heatMapZscoreSCRAMalgo(zscored_all_varcr, normiceall, bstring, varstring, stringfile, edges, limi, g)` | Builds and orders response heatmaps of z-scored neuronal activity around each stimulus/withdrawal event, returning the sorting order (`order`) of neurons used for visualization. |
| `cellClassBootStrapRAMalgo(varstring, bstring, normice, mice, files, edges)` | Classifies each neuron as activated or inhibited for each event type using a bootstrap-based statistical test, returning `actcells_all` and `inhcells_all` logical/index sets. |
| `heatmapResponsivenessMatriceOfEvents2(zscored_all_varcr, varstring, bstring, stringfile, actcells_allcr, inhcells_allcr, g)` | Builds a combined responsiveness matrix/heatmap across all event types, returning the response panel, combined percentages, and left/right response indices. |
| `shadedErrorBar(x, y, err, 'lineProps', ..., ...)` | Third-party plotting utility (Rob Campbell, MATLAB File Exchange) that draws a mean trace with a shaded error band (here, mean ± SEM) around it. |
| `smooth(y, span, method)` | Built-in MATLAB (Curve Fitting Toolbox) function used to smooth the mean and SEM traces before plotting, using a **loess** (local regression) method with span `smoothParam`. |

## Figures produced

- **Figure 1 — Z-scored activity**: grid of subplots with cell categories (`sensory`, `UR`, `PNLA`, `MNLA`, `TLA`, `MLA`, `PLA`, ...) as columns and selected events (thermal withdrawal, mechanical withdrawal) as rows. Each subplot overlays ipsilateral (paw stimulation ipsi to the recorded dorsal horn) vs. contralateral mean z-scored activity ± SEM.
- **Figure 2 — Normalized activity**: same layout, using normalized (non z-scored) firing rate instead of z-score.

