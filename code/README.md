# ELE419 DoS Laboratory — MATLAB project

Simulation and analysis of **Denial-of-Service (packet-drop) attacks** on a
scalar TCP-like networked LQG control loop, for the ELE419 lab brief
(`ELE419_DoS_lab_26.pdf`). Implements both exercises end-to-end with Monte Carlo
support, reproducible seeding, and automatic figure/table export.

## How to run

From a clean MATLAB session:

```matlab
cd '...\ELE419_DoS_matlab'      % the folder containing main_ELE419_DoS_lab.m
main_ELE419_DoS_lab
```

That single command creates `outputs/`, runs Exercise 1 and Exercise 2, saves
every figure (PNG + FIG) and table (CSV + MAT + LaTeX), writes a timestamped run
log, and prints the conclusions to the command window. Runtime ≈ 40 s on a
laptop (R2025b).

Tested on **MATLAB R2025b**. **No toolboxes required** — the simulation core
uses only base MATLAB (`randn`, `rand`). The shared template used `ss`,
`mvnrnd`, `binornd` (Control System + Statistics toolboxes); those are replaced
here by toolbox-free equivalents (`rand < p` for Bernoulli, `sqrt(var)*randn`
for Gaussians, plain structs instead of `ss`).

### Configuration (top of `main_ELE419_DoS_lab.m`)

```matlab
cfg.masterSeed = 2026;       % all randomness derives from this
cfg.mcTrials   = 300;        % Exercise 1 Monte Carlo trials per case
cfg.ex2grid    = 0:0.05:1;   % Exercise 2 (nu_bar, gamma_bar) grid (21x21)
cfg.ex2trials  = 80;         % Exercise 2 MC trials per grid point
cfg.ex2horizon = 100;        % Exercise 2 horizon (brief N=100)
```

Use `ex2grid = 0:0.02:1` for a finer (slower) map; lower `mcTrials` for a quick
dry run.

## Files

| File | Role |
|------|------|
| `main_ELE419_DoS_lab.m` | One-command driver: folders, seeding, both exercises, logging. |
| `run_exercise1.m` | Six DoS cases (A=1.1) with Monte Carlo; per-case figures + summary table. |
| `run_exercise2.m` | Stability search (A=1.255): deterministic recursions + Monte Carlo grid/sweeps. |
| `simulate_tcp_case.m` | Simulates **one** realisation of the scalar TCP-like LQG loop. |
| `kalman_tcp_update.m` | Clean, toolbox-free TCP-like Kalman update (replaces `KalmanTCP.m`). |
| `tcp_control_gain.m` | Steady-state TCP-LQG control gain from the modified Riccati (MARE). |
| `tcp_expected_covariance.m` | Expected error-covariance recursion (estimator MS-stability). |
| `stability_flag.m` | Heuristic finite-horizon stability label for Exercise 1. |
| `base_params.m` | Common configuration struct (notation matches the brief). |
| `plot_case_results.m` | Standard four-panel per-case diagnostic figure. |
| `make_summary_tables.m` | Exports a struct array to CSV + MAT + LaTeX. |
| `save_figure.m` | Robust PNG/FIG saver; forces a light (print-friendly) theme. |
| `README.md` | This file. |
| `report_support.md` | Report scaffold with the real generated numbers, captions, discussion. |
| `appendix_code_listing.md` | How to assemble the code appendix. |

## Outputs

```
outputs/
  exercise1/   case1..case6  (.png + .fig)          per-case 4-panel figures
  exercise2/   ex2_stability_heatmaps, ex2_cost_covariance_heatmaps,
               ex2_threshold_sweeps, ex2_deterministic_recursions (.png+.fig),
               exercise2_grid.mat
  figures/     ex1_cost_overlay, ex1_rms_state_bar  (cross-case comparisons)
  tables/      exercise1_summary.{csv,mat,tex}, exercise2_thresholds.{csv,mat,tex},
               all_results.mat
  run_log_YYYYMMDD_HHMMSS.txt
```

## Notation (matches the brief)

`nu_bar` = actuation arrival prob `ν̄`; `gamma_bar` = sensing arrival prob `γ̄`;
`Q` = process-noise variance `Var(w)`; `R` = measurement-noise variance
`Var(z)`; `Qxx, Quu` = LQ cost penalties (= 1). Time `k = 0..N` maps to MATLAB
index `k+1` (states `x_0..x_N`; inputs/measurements `k = 0..N-1`).

---

## Validation checklist (verify after you run)

- [ ] Exercise 1 uses `A=1.1, N=100, Var(w)=2`; the six cases match the brief's
      `(ν̄, γ̄, Var(z))`. → check the header line and Table 1.
- [ ] Exercise 2 uses `A=1.255, N=100, Var(w)=0.5, Var(z)=0.15`.
- [ ] Realised mean arrivals match the targets (e.g. Case 3 ≈ 0.05, Case 5 = 0).
- [ ] Case 1 has **zero** RMS estimation error (R=0 ⇒ exact estimate).
- [ ] Cases 3 & 5 keep small RMS estimation error (sensing intact); Cases 4 & 6
      have large RMS estimation error (sensing denied).
- [ ] Exercise 2 deterministic thresholds ≈ **0.365** = `1 − 1/A²` (Table 2).
- [ ] Re-running reproduces identical numbers (fixed seed).
- [ ] Figures saved under `outputs/`; tables under `outputs/tables/`.
- [ ] Original template values (`N=250, A=1.4, Q=0.5, R=0.1, nu_bar=gamma_bar=0.95`)
      do **not** appear in the assignment results.

## Three-lead review checklist

**Lead 1 — Control theory & assignment compliance.**
Checked: model (1)–(2), cost (3) with the `ν_k`-gated input penalty, TCP
information set `F_k={y^k,γ^k,ν^{k-1}}`, TCP Kalman prediction (uses known
`ν_{k-1}`) and γ-gated correction, steady-state TCP-LQG gain from the MARE, and
the `1−1/A²` threshold. Both exercises answered with the exact brief parameters.
Conclusions distinguish actuation loss from sensing loss (Cases 3↔4, 5↔6) and
explain the TCP separation that makes the Ex-2 stable region a rectangle.
*Assumption:* steady-state (not time-varying finite-horizon) gain — justified as
equivalent over `N=100`; **state this in the report.**

**Lead 2 — Simulation & reproducibility.**
Checked: one command reproduces all outputs; single master seed with
deterministic per-case/per-trial sub-seeds; Monte Carlo means ± std; base-MATLAB
only (toolbox-free, with `ss/mvnrnd/binornd` replaced); `inv()` avoided
(`/`, `\`); `R=0` handled safely; figures/tables auto-saved; light print theme.
*You must verify:* the project runs on **your** MATLAB (run it once); confirm
runtime is acceptable; optionally rerun Ex-2 with a finer grid.

**Lead 3 — Evidence & communication.**
Checked: every number in `report_support.md` traces to `outputs/tables/`;
discussion is specific ("realised ν̄ = 0.051", "RMS error ratio ≈ 0.31 vs 7×10³")
not generic; limitations are explicit (finite horizon, heavy-tailed second
moment, scalar-specific threshold); captions provided.
*You must do:* rewrite the discussion **in your own words**, trim to 7 pages,
choose which figures to show, and confirm the references.

### Remaining assumptions / things to confirm
- Initial state `x_0 ~ N(0, 0.01)` and prior mean 0 (matches the template's
  `Po=0.01`); change `p.x0`/`p.P0` in `base_params.m` if your reading differs.
- Stability flags in Ex-1 are heuristic; the rigorous boundary is Ex-2.
- The `1−1/A²` identity is exact for this scalar plant; cite Sinopoli (2004).
