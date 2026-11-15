# Error definitions and stored histories

The code preserves the manuscript's current definitions. It does not apply the alternative healthy-stiffness renormalization discussed during drafting.

## Final-ten-second stiffness NRMSE

For story i and case r, over the final time window W:

`100 * sqrt(mean(((estimated_k_i(t) - true_k_i(t)) / true_k_i(t))^2))`.

True stiffness is constant and post-damage in this window. The case metric averages the story-wise NRMSEs; the table averages those case metrics over successful runs. This is not endpoint MAE and not a pooled RMS across all stories. The final window is 93.740–103.740 s, inclusive: **10,001 samples** at 1 ms spacing.

## Response NRMSE and absolute spatial RMSE

For each displacement, velocity or acceleration channel in a window:

`NRMSE (%) = 100 * sqrt(sum((estimate - truth).^2) / sum(truth.^2))`.

The code averages channel NRMSEs within a case, then averages successful cases. A zero true-signal energy produces NaN. Very small true responses can give large relative errors even when physical errors are small.

`CASE_METRICS.csv` also reports the absolute spatial RMSE:

`sqrt(sum_over_channels_and_times((estimate - truth).^2) / (number_of_channels * number_of_times))`.

Its units are m, m/s or m/s^2. These quantities are distinct from the channel-mean NRMSEs.

Windows retain the archived sample indices: full record; healthy forced response after 5 s; first 5 s after the damage switch; remaining forced response; free decay; and final 10 s. All metric accumulators use every update, not only the saved plotting samples.

## KDE + box distributions

Stiffness uses signed healthy-normalized error:

`100 * (estimated_k_i(t) - true_k_i(t)) / healthy_k_i`.

Displacement, velocity and acceleration use signed **absolute** errors `estimate - truth`, in m, m/s and m/s^2. These are not NRMSEs. Each story pools 94 integer-second samples (10–103 s) per successful proposed-method case: 9,400 samples for a complete 100-case ensemble.

The Gaussian KDE sums normal kernels centered at every observation, divided by sample count and bandwidth. Width is rescaled graphically to equal peak width per story. Boxes contain the 25th–75th percentiles; the median is marked. Whiskers end at the outermost actual samples within 1.5 IQR of the box boundaries. They are not confidence limits. Central ranges retain the whiskers; full-range companion PDFs include all observations. Temporal samples are correlated and are not independent Monte Carlo replications.

To avoid drawing a horizontal baseline across every story, the visible violin outline ends where density falls below 0.1% of that story's peak. This affects only the drawn outline; all observations still enter the KDE, quantiles and extrema, and the full-range axes cover all observations.

## NIS and timing

The single-case plot shows empirical NIS divided by measurement count (39 for Pirelli, 4 for ASCE), with a reference line at one. This empirical statistic uses the learned innovation baseline. It is not itself a calibrated damage alarm. The applied forgetting factor is captured at every 1 ms update.

The trace contains 0.1 s samples plus boundary/end samples for stiffness, conditional posterior stiffness SD, response estimates and truth. Full-resolution arrays are retained only for applied alpha, empirical normalized NIS, severity and core update time. Large complete response histories, observations and process-covariance histories are not stored.

Single-case stiffness plots show conditional posterior +/-2 SD bands. These are uncertainty estimates from the filter, not independently validated confidence or coverage intervals.
