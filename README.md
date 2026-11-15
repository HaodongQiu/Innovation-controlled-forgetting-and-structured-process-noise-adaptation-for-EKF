# Adaptive EKF for joint structural response and stiffness estimation

MATLAB code and compact input data accompanying **Innovation-controlled forgetting and structured process-noise adaptation for joint state–stiffness estimation in structural dynamics**, by Haodong Qiu, Aldo Ghisi and Stefano Mariani.

The release contains the six configurations in the manuscript, the 39-DOF equivalent Pirelli Tower model, the four-story IASC–ASCE shear-building model, the supplied El Centro earthquake input, and the matched 100-case experimental designs. Simulation and plotting are separate.

## Requirements

- The paper campaigns used MATLAB **R2024b Update 1** on Windows. Local release validation is recorded in `docs/VALIDATION.md` with the actual version used.
- Base MATLAB supports serial calculations and all plotting, including the Gaussian KDE and box statistics. No Statistics and Machine Learning Toolbox is needed.
- Parallel Computing Toolbox is required only for runs with more than one process worker. Four workers are recommended; use `1` for serial execution.
- CPU arrays are used. A GPU and CUDA are not required.

## Start with one filter

Download/clone the repository, set MATLAB's Current Folder to the repository root, and run:

```matlab
setup_project
result = run_single_filter('Pirelli39','M12-dynamic',51);
plot_single_result(result)
```

This uses case 51, measurement-noise SD `0.01 m/s^2`, and Newmark-generated responses. It saves one result under `outputs/single/` and generates vector PDFs only when `plot_single_result` is called. Case 51 is a favorable illustration selected by the smallest final stiffness NRMSE; it does not represent typical ensemble accuracy.

For the four-story benchmark or another configuration:

```matlab
result = run_single_filter('ASCE4','U-fixed',1);
result = run_single_filter('Pirelli39','M1-dynamic',1);
result = run_single_filter('Pirelli39','M12-dynamic',51,0.03,'ode45');
```

## Six configurations

Use these exact method IDs in commands. `paper_methods` is the single configuration catalogue.

| ID | Response covariance mapping (M1) | Independent parameter disturbances / zero response–parameter process covariance (M2) | Response forgetting | Stiffness forgetting |
|---|---|---|---|---|
| `U-fixed` | No | No | 0.80 | 0.80 |
| `U-dynamic` | No | No | 0.95 | Dynamic, common to all parameters |
| `M1-dynamic` | Yes | No | 0.95 | Dynamic, common to all parameters |
| `M2-dynamic` | No | Yes | 0.95 | Dynamic, common to all parameters |
| `M12-fixed` | Yes | Yes | 0.80 | 0.80 |
| `M12-dynamic` | Yes | Yes | 0.95 | Dynamic, common to all parameters; proposed |

“Fixed” refers to the forgetting factor: process-noise covariance adaptation remains active in all six configurations. The dynamic stiffness factor starts at 0.80 and is bounded by 0.80 and 0.95. There is no localization. M2 restricts the **process-noise** covariance; the posterior state covariance remains full. Fixed/dynamic comparisons also differ in the response forgetting factor, as in the manuscript; they do not isolate only the stiffness schedule.

## Matched 100-case Monte Carlo studies

```matlab
run_monte_carlo(4)  % both buildings, six methods, 100 cases each: 1,200 outcomes
```

Each case generates truth once, then applies every requested method to identical truth and identically seeded measurement noise. The Pirelli design has three damaged stories per case; ASCE has one. Reductions are between 1% and 20%. Frozen case files reproduce the paper's exact designs instead of drawing a new ensemble.

Useful alternatives:

```matlab
run_monte_carlo(1)  % serial; no parallel toolbox
run_monte_carlo(4,{'Pirelli39'})
run_monte_carlo(4,{'ASCE4'},0.01,'newmark',1:100,{'M12-dynamic'})
run_monte_carlo(1,{'ASCE4'},0.01,'newmark',1:2,{'M12-dynamic'}) % small trial
```

Completed checkpoints are reused when the **same command and unchanged code/data** are used again. A failed filter is saved as a completed failed comparison and is not repeatedly rerun. A SHA-256 signature prevents reuse after source/data changes; move old outputs aside if you deliberately change the experiment. Temporary files are excluded from summaries. Avoid concurrent runs writing to the same campaign directory.

## Noise and response-generation studies

```matlab
run_noise_sensitivity(4)        % proposed Pirelli filter: SD 0.005, 0.01, 0.03, 0.05 m/s^2
run_response_generator_study(4) % proposed Pirelli filter: Newmark, ode45, matrix exponential
```

The filter always uses its original explicit Newmark transition. The alternative response generators change only the simulated truth. `ode45` and `expm_foh` use the same piecewise-linear input; their implementation and tolerances are in `src/continuous_response.m` and `src/paper_settings.m`.

To run all manuscript campaigns:

```matlab
run_paper_studies(4)
```

This produces **1,700 unique outcomes**: 1,200 benchmark comparisons, 300 additional proposed noise cases and 200 additional response-generator cases. The reference-noise Newmark cases are reused. This is a substantial calculation, not a quick example.

## Results and plotting

The simulation commands save results and numerical summaries; they do not create figures.

```matlab
plot_paper_results  % generate vector PDFs from all saved Monte Carlo campaigns
```

Or select one saved campaign:

```matlab
folder = campaign_path('Pirelli39',0.01,'newmark');
summarize_results(folder) % regenerate CSV, JSON and Markdown only
plot_paper_results(folder)
```

Each campaign contains finalized `case_####_METHOD.mat` files, `CASE_METRICS.csv`, `METHOD_SUMMARY.csv`, `summary.json`, `RESULTS_REPORT.md`, and `pipeline_status.json`. PDFs are written to its `figures/` subfolder. `summary.complete` means 100 finalized outcomes for every method present in that campaign; a requested small trial remains a partial ensemble.

Plots include comparison CDFs, endpoint estimated/true stiffness ratios, story-wise final-ten-second stiffness NRMSE, and horizontal KDE + box distributions for stiffness, displacement, velocity and acceleration. The single-case plotter also exports empirical NIS normalized by the measurement count and the applied forgetting factor. Plotting code is exclusively under `plots/`, apart from the root plotting entry point.

KDE plots use 1 s samples at 10, 11, ..., 103 s (94 times per successful case). The filter updates and error accumulators use every 1 ms sample. Each plot has a central-view PDF retaining all Tukey whiskers and a full-range companion. All samples enter the KDE and statistics. Numeric extrema and bandwidths are saved in `KDE_METHOD.json`; samples are correlated in time. For Pirelli the bandwidths are 0.03 percentage points for stiffness, 7.3e-6 m, 2.1e-5 m/s and 3.6e-4 m/s^2 for the three responses.

The plotters regenerate the scientific figure types from new results; they do not reproduce the manuscript's complete LaTeX layout or every cosmetic annotation. See [metric definitions](docs/METRICS.md) before comparing stiffness and response distributions.

## Data and settings

See [data sources and variable descriptions](data/README.md), [numerical settings](docs/SETTINGS.md), [references](docs/REFERENCES.bib), and [source provenance](docs/SOURCE_PROVENANCE.json).

The original archived sampling/index convention is preserved: `dt=0.001 s`, 53.740 s of earthquake input, and a 103.740 s record including free decay. The first damaged sample is at 26.870 s; the manuscript displays damage onset as 27 s. No archived data or case timing was shifted for this release.

The 20-DOF exploratory study, obsolete localized controllers, fixed-0.95/0.50 supplements, dynamic lower-bound 0.10 experiment, 500-case campaigns, bulky response histories, old output folders and manuscript drafts are intentionally outside this release.

## Validation and citation

```matlab
setup_project
validate_release
```

This checks input files and syntax and runs 12 short model/configuration checks. It does not launch a full Monte Carlo campaign. Release-specific regression evidence is in `docs/VALIDATION.md`.

The rounded manuscript results at the reference noise level are recorded in `docs/PAPER_REFERENCE_RESULTS.csv`. Accuracy averages exclude failed runs; numerical failure times can vary slightly with MATLAB/BLAS versions. These reference numbers are supplied for checking, not substituted into newly computed summaries.

Please cite the accompanying manuscript and the original structural-data references when using this code/data. `CITATION.cff` contains the manuscript metadata. Publication DOI and a software license have not been supplied; no license for third-party source data is inferred from their inclusion.
