# Input data

All numerical simulation arrays use SI units: mass kg, stiffness N/m, time s, displacement m, velocity m/s and acceleration m/s^2. These are model inputs and a supplied earthquake record, not new experimental measurements collected for this release.

| File | Variables | Purpose |
|---|---|---|
| `ElCentro_time_history_dt0001.mat` | `time_new`, `acc_new` (53,741 samples) | Supplied El Centro acceleration, at 0.001 s spacing. No g-to-SI conversion or amplitude rescaling is performed by this release. |
| `models/Pirelli39_model.mat` | `M`, `K0`, `k0`, `response_acc_variance`, `n`, `N`, `half_N`, `forced_samples`, `dt`, `zeta` | Simulation-ready equivalent 39-DOF shear-building model and prescribed priors. |
| `models/ASCE4_model.mat` | Same fields, four DOFs | Four-story shear-building representation of the IASC–ASCE steel-frame benchmark. |
| `cases/Pirelli39_cases.mat` | `case_id`, `true_k`, `damaged_stories`, `damage_pct`, `case_seed`, `noise_seed` | Frozen 100-case Pirelli design: three damaged stories per realization. |
| `cases/ASCE4_cases.mat` | Same fields | Frozen 100-case ASCE design: one damaged story per realization. |
| `reference/Pirelli_reference_mass_stiffness.mat` | `My`, `Ky` | Supplied original reference matrices (`KyMy.mat`), retained for traceability. |
| `reference/Pirelli_equivalent_story_stiffness.mat` | `K` (39-vector) | Supplied equivalent story stiffnesses (`K_final_Mode1.mat`). `K` here is a vector, not the assembled matrix. |

The filter uses the simulation-ready `M`, `K0` and `k0`, **not the original full reference stiffness matrix `Ky`**. With story-difference operator E, `K0 = E' * diag(k0) * E`. The Pirelli equivalent model was constructed to retain the reference first mode and fundamental frequency; equivalence of all higher modes is not asserted. The first three model levels are underground; levels 4–6 connect to the adjacent building.

The four-story reference benchmark is a steel frame; this study uses its four-story shear-building representation. See Johnson et al. (2004) for benchmark details. The exact matrices supplied to these simulations are retained in the model MAT file.

Pirelli reference data: Ebrahimzadeh Hassanabadi et al. (2020), DOI `10.1016/j.cma.2020.113334`; see also Qiu et al. (2024), DOI `10.3390/ecsa-11-20493`. ASCE reference: Johnson et al. (2004), DOI `10.1061/(ASCE)0733-9399(2004)130:1(3)`. Full bibliographic entries are in `docs/REFERENCES.bib`.

The supplied El Centro MAT file has no embedded station/component, preprocessing history or redistribution-license metadata. It is included unchanged as the exact numerical input used for these simulations. Do not infer additional acquisition details from its filename.

Case columns correspond to case IDs 1–100. Story numbering is one-based from the bottom. `true_k` is the full post-damage stiffness vector; unchanged entries equal `k0`. `damage_pct` contains the imposed reductions in percent, not fractions. `case_seed` records the original design provenance; simulation uses the stored designs directly. `noise_seed` drives a local `mt19937ar` stream, drawing one full channel at a time, identically across methods and across noise levels.

No simulated full-response history is needed as an input. `newmark_response` or `continuous_response` regenerates it from the compact model, earthquake and case design. Redundant initial covariance arrays and saved reference responses have been removed from the public case files; the same initial quantities are computed by the filter.
