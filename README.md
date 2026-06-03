# DoS Attacks on Networked Control Systems

Simulation and stability analysis of **Denial-of-Service (packet-drop) attacks** on a
scalar TCP-like **networked LQG control loop**. The project models an unstable plant
controlled over a lossy network and quantifies how dropping packets on the **sensing**
channel versus the **actuation** channel degrades estimation and control — including the
exact mean-square stability boundary for the closed loop.

Implemented end-to-end in **base MATLAB** (no toolboxes), with reproducible seeding and
automatic figure/table export. Accompanied by a written report (PDF) presenting the
results.

> Coursework for ELE419 *Cybersecurity for Control Systems*. All code and the report in
> this repository are my own work.

## Problem

A scalar discrete-time plant `x_{k+1} = A x_k + ν_k u_k + w_k` is observed and actuated
over a network where packets arrive with probabilities `γ̄` (sensing) and `ν̄`
(actuation). The controller/estimator is the TCP-like LQG (Kalman estimator + modified
Riccati / MARE control gain), where dropped packets are acknowledged. Two questions are
studied:

- **Exercise 1 (`A = 1.1`, `N = 100`, `Var(w) = 2`).** Six DoS cases compared. The key
  contrast is *which channel is attacked*: losing **actuation** keeps the state estimate
  accurate but loses control authority, while losing **sensing** blows the estimate up.
- **Exercise 2 (`A = 1.255`, `N = 100`, `Var(w) = 0.5`, `Var(z) = 0.15`).** Find the
  region of arrival probabilities for which the loop is mean-square stable. Result: the
  loop is **MS-stable iff `ν̄ > 0.365` and `γ̄ > 0.365`** — a rectangle, reflecting the
  TCP separation principle. The threshold matches the analytic critical probability
  `1 − 1/A² = 0.365`.

A subtle point captured in the analysis: Monte-Carlo sample means *under-estimate* the
mean-square boundary (the second moment is heavy-tailed), so the precise threshold is
obtained from **deterministic recursions** (control MARE + expected-covariance
recursion). The Monte-Carlo grid is illustrative.

## Repository layout

```
.
├── code/        MATLAB source (run main_ELE419_DoS_lab.m) + figures + notes
│   ├── *.m
│   ├── figures/         exported result figures (PNG)
│   └── README.md        detailed code documentation, parameters, validation checklist
├── report/      written report
│   ├── Nijat_Badalov_ELE419_DoS_Lab_Report.pdf
│   ├── Appendix_Code_Exercise1.m   self-contained reproduction of Exercise 1
│   ├── Appendix_Code_Exercise2.m   self-contained reproduction of Exercise 2
│   └── figures/         figures used in the report (PNG)
└── README.md
```

## Running

From a clean MATLAB session (tested on **R2025b**, no toolboxes required):

```matlab
cd code
main_ELE419_DoS_lab
```

One command creates `outputs/`, runs both exercises, saves every figure and table,
writes a timestamped run log, and prints the conclusions. Runtime ≈ 40 s on a laptop.
All randomness derives from a single master seed (`2026`), so re-running reproduces
identical numbers. See [`code/README.md`](code/README.md) for parameters, file roles,
and a validation checklist.

The two `report/Appendix_Code_Exercise*.m` scripts are standalone — each reproduces the
numbers for one exercise on its own.

## Key results

| | Exercise 1 (`A = 1.1`) | Exercise 2 (`A = 1.255`) |
|---|---|---|
| Critical arrival prob `1 − 1/A²` | 0.174 | 0.365 |
| Actuation loss (sensing intact) | small RMS estimation error, control lost | — |
| Sensing loss | estimate diverges | — |
| MS-stable region | — | `ν̄ > 0.365` **and** `γ̄ > 0.365` |

## License

Code is released under the [MIT License](LICENSE). The report PDF is academic work,
included for reference.
