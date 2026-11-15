# Numerical settings

`src/paper_methods.m` defines the six method IDs. `src/paper_settings.m` supplies response factors, metric windows and solver tolerances. `src/aekf_filter.m` retains the covariance update and numerical safeguards of the paper implementation.

| Setting | Value |
|---|---|
| Time step | 0.001 s |
| Newmark parameters in filter and default truth | beta = 0, gamma = 0.5 |
| Rayleigh damping | Fixed from healthy model; damping ratio 0.01 at the first two modes |
| Initial response estimate | Zero |
| Initial stiffness estimate / healthy stiffness | 1.50 |
| Initial stiffness SD / healthy stiffness | 0.50 |
| Prescribed acceleration-prior variances | `response_acc_variance` in each model file |
| Initial velocity-prior variance | max((10 dt)^2 v_a, 1e-6) m^2/s^2 |
| Initial displacement-prior variance | max((0.5 (10 dt)^2)^2 v_a, 1e-10) m^2 |
| Initial acceleration driving variance | 10 v_a |
| Acceleration driving-variance floor | 0.01 times the measurement variance |
| Stiffness diffusion amplitudes | Lower 0.001, initial/nominal 0.01, upper 0.05 s^(-1/2) |
| Stiffness driving variances | (amplitude times healthy stiffness)^2 dt |
| Fast / baseline / warm-up innovation factors | 0.98 / 0.9995 / 0.99 |
| Warm-up duration | 5 s |
| Change threshold / ramp width | 1.15 / 1.50 |
| Dynamic stiffness-factor bounds | 0.80 to 0.95, initially 0.80 |
| Attack / recovery rates | 0.20 / 0.002 per sample |
| Fixed configurations | Both response and stiffness factors equal 0.80 |
| Dynamic configurations | Response factor fixed at 0.95 |
| Gaussian measurement-noise SD | Default 0.01 m/s^2; noise study uses 0.005, 0.01, 0.03, 0.05 |
| Measurement covariance | Diagonal, using each generated noise channel's sample variance, as in the archived study |
| ode45 tolerances | RelTol 1e-8; displacement AbsTol 1e-11 m; velocity AbsTol 1e-10 m/s; MaxStep 0.001 s |

These prescribed initialization quantities do not use the current realization's future true response. Measurement covariance is computed from the generated noise in this synthetic experiment; an application to actual measurements would need its own sensor-noise calibration.

The filter uses one linear-algebra thread per MATLAB process and restores the caller's thread setting afterwards. Parallelism is across independent cases, not across time updates. `core_step_seconds` records update time excluding error accumulation and passive history capture; it includes timing-call overhead. The observed mean update time is not a hard real-time latency guarantee.

No stiffness projection to positive values is added. A run is successful when the entire record is processed without divergence, all final states are finite and all final stiffnesses are positive. Covariance-factorization failures, overflow/nonfinite states and nonpositive final stiffness produce failed outcomes. No arbitrary accuracy threshold is used as a success criterion.
