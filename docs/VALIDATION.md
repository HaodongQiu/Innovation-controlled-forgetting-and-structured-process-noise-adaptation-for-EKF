# Release validation

The publication campaigns used MATLAB **R2024b Update 1**. The local installation available for release checks is **R2023b Update 7 (23.2.0.2515942)**, using CPU arrays and one numerical thread per filter. A complete 100-case rerun was not performed for packaging.

## Source and input checks

- The public filter was derived from the supplied `independent_filter.m`, which supports the same six covariance/forgetting configurations. After removing the new passive history/timing assignments, its active computational statements recover the original source exactly. Equations, covariance projections, numerical safeguards and random-noise generation were not retuned.
- Pirelli `M` is identical to the supplied `My`, and `k0` is identical to the supplied equivalent story-stiffness vector. Both models satisfy the supplied shear-matrix assembly relation.
- The compact case files retain the original 100 damage designs, one-based story IDs, damage magnitudes and seeds. Only redundant reference responses and repeated initialization arrays were removed.
- All MATLAB function names match their files. Simulation code contains no plot-generation calls or absolute workstation paths.

## Executed local checks

1. **Twelve short smoke tests:** all six configurations on both structures, 600 updates each. Input dimensions, finite output, fixed alpha = 0.80 and dynamic alpha bounds were checked.
2. **Six same-platform source comparisons:** original and cleaned filters were run on the same four-story inputs for 6,200 updates, including a damage transition after controller warm-up. Every configuration produced **zero difference** in final states and accumulated response-error energies. Details: `SAME_PLATFORM_REGRESSION.json`.
3. **Alternative response generators:** on a shortened four-story case, maximum `ode45`/matrix-exponential differences were 1.30e-12 m, 3.17e-10 m/s and 6.79e-8 m/s^2. Details: `PARALLEL_AND_TRUTH_VALIDATION.json`.
4. **Checkpoint reuse:** two complete four-story proposed-method cases were run at a diagnostic noise SD of 0.0123 m/s^2; repeating the command reused both finalized results. This diagnostic noise is outside the paper's four-level study.
5. **Full-record pilots:** case 1 was processed with all six configurations for both structures (12 finalized, successful pilot outcomes). This is not a replacement for the paper's 100-case ensembles.

6. **Plot export and visual review:** 22 pilot PDFs (Pirelli case 1) were exported, rendered with Poppler and inspected. The KDE plots in this smoke test pool 94 observations per story from one case; production plots pool 9,400 after all 100 successful cases are available. CDFs also display correctly for one-case trials. See `PLOT_VALIDATION.json`.

Parallel Computing Toolbox is **not installed locally**. Parallel execution was therefore not executed here. The process-pool runner is provided for a workstation with the toolbox; installation and license are checked before creating a pool. Use `workers=1` for serial operation.

## Cross-version reproducibility

Full-record pilot outcomes are compared with the stored R2024b paper archives in `REGRESSION_RESULTS.json`. An `archive_close` flag records agreement within the comparison tolerance; it must not be read as universally true. The correlated-covariance configurations exhibit differences between the local pilots and those archives despite identical initialization and damage data and measurement variances agreeing to floating-point precision.

The same-platform source tests above distinguish a cleanup error from cross-platform/version reproducibility. They establish equality with the original implementation on the tested local inputs. They do not prove bitwise equality with every MATLAB/BLAS implementation or identical failure counts over a new full ensemble. Use the paper's MATLAB release for the closest reproduction. The supplied `PAPER_REFERENCE_RESULTS.csv` retains the manuscript's results and is never used to replace computed outputs.

For the full-record Pirelli case-1 pilots, M2-dynamic, M12-fixed and M12-dynamic have maximum scaled final-state differences of approximately 1.28e-13, 1.76e-13 and 9.65e-14 against their archives. U-fixed, U-dynamic and M1-dynamic differ more, despite identical initial states and covariances; their individual comparisons remain visible in the JSON report rather than being suppressed.

All six ASCE case-1 final-state comparisons agree within 4.3e-12 in the scaled norm, and their story-wise final stiffness NRMSEs agree within 4.6e-10 percentage points. The final response NRMSE comparisons are more sensitive because the true free-decay response is very small. Accordingly, the strict combined `archive_close` flag can be false even when the ASCE stiffness and final-state comparisons agree closely.

Generated validation cases and PDFs are kept outside the public source package. The `outputs/` directory is created when a user runs the examples and is excluded by `.gitignore`.
