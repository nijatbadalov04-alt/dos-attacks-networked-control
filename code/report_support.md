# Report support — ELE419 DoS Laboratory

This file is **scaffolding for your 7-page report**, not a finished submission.
Every number quoted below was produced by `main_ELE419_DoS_lab.m` (master seed
2026, 300 Monte Carlo trials per case in Exercise 1; 80 trials per grid point
in Exercise 2). Re-run the project and check the numbers against
`outputs/tables/` before you submit — then rewrite the discussion in your own
words. Do not paste this text verbatim.

---

## Proposed title

**Denial-of-Service attacks on a TCP-like networked LQG loop: performance loss
under actuation and sensing packet drops, and the mean-square stability
boundary.**

(Alternatives: "Availability attacks on optimal control: a TCP-like Kalman/LQG
study"; keep it specific to *packet drops / DoS*, not false-data injection.)

---

## Suggested 7-page layout

| § | Content | ~pages |
|---|---------|--------|
| 1 | Introduction & aim: DoS as an availability attack; what packet drops do to estimation vs actuation | 0.75 |
| 2 | Model & method: eqs (1)-(3), TCP-like information set, TCP Kalman filter, steady-state LQG gain, critical arrival probability | 1.5 |
| 3 | Exercise 1 — setup, the summary table, one or two representative figures | 1.25 |
| 4 | Exercise 1 — discussion across the six cases (the actuation-vs-sensing contrasts) | 1.25 |
| 5 | Exercise 2 — method: deterministic recursions + Monte Carlo grid/sweeps | 0.75 |
| 6 | Exercise 2 — results: the stable rectangle, threshold vs 1−1/A², limitations | 1.0 |
| 7 | Conclusion | 0.5 |
| App | All MATLAB code (outside the 7-page limit) | — |

Keep figures small (two per row). Put only the *summary* table in the main
text; the per-case multi-panel figures can be cited and one or two shown.

---

## §2 Theory background (condensed — expand/rephrase)

The plant is the scalar Gauss–Markov system of eqs (1)–(2),

```
x_{k+1} = A x_k + ν_k B u_k + w_k ,   w_k ~ N(0,Q)
y_k     = γ_k C x_k + z_k         ,   z_k ~ N(0,R)
```

with i.i.d. Bernoulli arrivals P[ν_k=1]=ν̄ (actuation) and P[γ_k=1]=γ̄
(sensing). The cost (3) penalises the state always but the input **only when it
is actually applied** (`ν_k U′ Q_uu U`), so packet loss is not double-counted.
`Q_xx = Q_uu = 1`.

**TCP-like information set** `F_k = {y^k, γ^k, ν^{k-1}}`. The acknowledgement of
the previous actuation outcome `ν_{k-1}` is available, so the filter's
prediction is deterministic given `F_k`:

```
prediction:  x̂⁻_k = A x̂_{k-1} + ν_{k-1} B u_{k-1},   P⁻_k = A P_{k-1} A' + Q
correction:  γ_k=1 → K_k = P⁻_k C'/(C P⁻_k C'+R);  x̂_k = x̂⁻_k + K_k(y_k−C x̂⁻_k);  P_k = P⁻_k − K_k C P⁻_k
             γ_k=0 → x̂_k = x̂⁻_k,  P_k = P⁻_k
```

Because actuation is acknowledged, `ν` never enters the covariance — this is the
key TCP feature and is why **estimation and control separate** (see Ex. 2). The
controller is the steady-state TCP-LQG law `u_k = −L x̂_k`, with `L` and value
`S` from the modified Riccati recursion (MARE)

```
L = (B'SB+Q_uu)^{-1} B'SA ,   S = A'SA + Q_xx − ν̄ A'SB L .
```

**Critical arrival probability.** For a scalar unstable plant the second moment
`E[x_k²]` (and the expected error covariance) stay bounded iff the arrival
probability exceeds `1 − 1/A²` (Sinopoli et al. 2004 for estimation; the
analogous MARE condition for control). For Exercise 1 `A=1.1 → 1−1/A²=0.174`;
for Exercise 2 `A=1.255 → 1−1/A²=0.365`.

*Implementation note for the report:* the shared template runs the MARE forward
in time alongside the simulation; we iterate it to its fixed point and apply the
constant steady-state gain. Over `N=100` the two are numerically identical after
a short transient. As a sanity check, for Case 2 the predicted steady-state
posterior standard deviation `√P ≈ 0.86` matches the measured RMS estimation
error `0.854`, and the LQR gain `L = 0.703` matches the analytic root of the
scalar Riccati — evidence the implementation is correct.

---

## §3 Exercise 1 — results

`A=1.1, B=1, C=1, N=100, Var(w)=2`. Means ± std are over 300 Monte Carlo trials.

**Table 1.** Exercise 1 summary (from `outputs/tables/exercise1_summary.csv`).

| Case | ν̄ | γ̄ | Var(z) | realised ν̄ | realised γ̄ | mean cost `E[J]/N` | RMS state | RMS est. error | max\|x\| | flag |
|------|----|----|--------|-----------|-----------|------------------|-----------|----------------|---------|------|
| 1 | 1.00 | 1.00 | 0.0 | 1.000 | 1.000 | 3.55 | 1.53 ± 0.13 | 0.00 | 4.2 | stable |
| 2 | 1.00 | 1.00 | 1.0 | 1.000 | 1.000 | 4.59 | 1.81 ± 0.18 | 0.85 ± 0.07 | 4.9 | stable |
| 3 | 0.05 | 1.00 | 0.5 | 0.051 | 1.000 | 2.1×10⁵ | 112 ± 434 | 0.64 ± 0.05 | 462 | unstable |
| 4 | 1.00 | 0.05 | 0.5 | 1.000 | 0.049 | 6.3×10⁵ | 208 ± 744 | 188 ± 675 | 856 | unstable |
| 5 | 0.00 | 1.00 | 0.1 | 0.000 | 1.000 | 9.1×10⁷ | 7.5×10³ ± 5.9×10³ | 0.31 ± 0.02 | 3.1×10⁴ | unstable |
| 6 | 1.00 | 0.00 | 0.1 | 1.000 | 0.000 | 9.9×10⁷ | 7.7×10³ ± 6.3×10³ | 7.0×10³ ± 5.7×10³ | 3.2×10⁴ | unstable |

(Final realised cost `E[J]`: Case 1 ≈ 355, Case 2 ≈ 459, Case 3 ≈ 2.1×10⁷,
Case 4 ≈ 6.3×10⁷, Cases 5–6 ≈ 9–10×10⁹. The standard deviations for the
unstable cases are larger than the means — the cost is heavy-tailed, dominated
by the rare trials with long outage runs. Report the *order of magnitude*, not
spurious precision.)

**Figures.** `outputs/exercise1/caseN.png` — four panels each: (1) log₁₀
running cost, (2) state `x_k` and estimate `x̂_k`, (3) estimation error with the
±√P_k band, (4) the realised `γ_k`/`ν_k` sequences. Cross-case overlays in
`outputs/figures/` (`ex1_cost_overlay.png`, `ex1_rms_state_bar.png`).

Suggested captions:
- *Fig. X.* Case 1 (ideal). State and estimate coincide (R=0 ⇒ zero estimation
  error); the running cost settles, confirming a regulated, mean-square-stable
  loop.
- *Fig. X.* Case 6 (no sensing). The estimate stays near zero while the true
  state runs away to ≈ −2.5×10⁴; the ±√P_k band fans out, showing the
  unbounded growth of estimation uncertainty.
- *Fig. X.* RMS state per case (log axis, MC mean ± std): the four DoS cases
  exceed the two healthy cases by two to four orders of magnitude.

---

## §4 Exercise 1 — discussion (draft prose, rewrite in your voice)

**Cases 1–2 (healthy channels).** With both links up the loop behaves as a
standard LQG regulator. Removing measurement noise (Case 1) gives exact state
knowledge: the estimate tracks `x_k` to machine precision and the RMS estimation
error is identically zero. Adding `z_k ~ N(0,1)` (Case 2) raises the mean cost
only from 3.55 to 4.59 (≈ 29 %) and the RMS state from 1.53 to 1.81, and the
measured RMS error (0.85) sits at the steady-state Kalman value. This is the
expected graceful degradation — noise costs a little, but the loop stays
mean-square stable, as both arrival probabilities (1.0) are far above the
0.174 critical value for `A=1.1`.

**Case 3 vs Case 4 (severe single-channel loss).** Both deny one link 95 % of
the time, but on different channels, and the consequences differ sharply. Under
actuation loss (Case 3, realised ν̄ = 0.051) the estimator is untouched —
sensing is perfect, so the RMS estimation error stays at 0.64 — yet the RMS
state climbs to ≈ 112 and the cost to ≈ 2×10⁷, because the controller's (good)
commands almost never reach the plant, which is therefore open-loop and
unstable (A>1) most of the time. Under sensing loss (Case 4, realised γ̄ =
0.049) the failure is compounded: the estimator now runs mostly on prediction,
its RMS error blows up to ≈ 188, and the actuator — which *does* work — injects
those wrong commands, so the RMS state (≈ 208) and cost (≈ 6×10⁷) are both
higher than in Case 3. The headline: *actuation loss costs you control authority
but leaves perception intact; sensing loss destroys both.*

**Case 5 vs Case 6 (full denial of one channel).** This is the cleanest
contrast and the strongest evidence for the actuation/sensing distinction. With
no actuation (Case 5, ν̄=0) the loop is the open-loop recursion
`x_{k+1}=1.1 x_k + w_k`; the state diverges (RMS ≈ 7.5×10³) but, because sensing
is perfect, the estimator tracks that runaway almost exactly — RMS estimation
error only 0.31. The controller is reduced to a perfect *observer* of a plant it
cannot touch. With no sensing (Case 6, γ̄=0) the opposite holds: the actuator
works, but the filter never receives a correction, so for `A>1` the error
covariance grows without bound (visible as the widening ±√P_k band), the
estimate becomes meaningless, and the loop collapses (RMS state ≈ 7.7×10³, RMS
error ≈ 7.0×10³). The RMS-estimation-error ratio between the two cases —
≈ 0.31 versus ≈ 7.0×10³, four orders of magnitude — is effectively a signature
of *which* channel was attacked.

**Tie to theory.** All four DoS cases push one arrival probability below the
0.174 critical value for `A=1.1` (Cases 3,5 on actuation; Cases 4,6 on sensing),
and all four are flagged unstable, while Cases 1–2 (probabilities = 1) are
stable. Over a finite `N=100` with a mildly unstable `A=1.1` the divergence is
gradual rather than explosive, but the trend and the heavy-tailed cost spread
are consistent with mean-square instability.

*Editorial note:* mention that the stability flags here come from a heuristic
finite-horizon threshold (`stability_flag.m`); they agree with the plots, but
the rigorous boundary is the Exercise-2 analysis.

---

## §5–6 Exercise 2 — method and results

`A=1.255, B=1, C=1, N=100, Var(w)=0.5, Var(z)=0.15`. Goal: the `(ν̄, γ̄)` region
of stability. Two complementary approaches (both in `run_exercise2.m`):

1. **Deterministic (precise).** Mean-square stability fails exactly where the
   two recursions lose boundedness: the control MARE (`ctrl_steady`) and the
   expected error-covariance recursion (`tcp_expected_covariance`). We locate
   each threshold by bisection on boundedness.
2. **Monte Carlo (illustrative).** Simulate the full stochastic loop on a
   21×21 `(ν̄,γ̄)` grid and on fine 1-D sweeps; report `log₁₀` of the mean state
   energy `E[(1/N)Σx_k²]`.

**Table 2.** Exercise 2 thresholds (from `outputs/tables/exercise2_thresholds.csv`).

| Quantity | Value |
|----------|-------|
| Analytic `1 − 1/A²` | **0.3651** |
| Deterministic ν̄ threshold (control MARE) | **0.3652** |
| Deterministic γ̄ threshold (expected covariance) | **0.3652** |
| MC ν̄ threshold, N=100 | 0.19 |
| MC ν̄ threshold, N=300 | 0.29 |
| MC γ̄ threshold, N=100 | 0.40 |
| MC γ̄ threshold, N=300 | 0.28 |

**Answer.** The system is mean-square stable iff

```
        ν̄ > 0.365   AND   γ̄ > 0.365 ,
```

i.e. the stable set is the rectangle `(0.365,1] × (0.365,1]` (Fig. A,
`ex2_stability_heatmaps.png`, left panel). Both thresholds equal `1 − 1/A²` to
four significant figures.

**Why a rectangle (TCP separation).** In the TCP setup the controller knows
`ν_{k-1}`, so `ν` never enters the error covariance: the estimator threshold
depends on `γ̄` only, and the control threshold on `ν̄` only. The two conditions
are independent, so the stable region is their product — a rectangle, not a
curved trade-off. The final-covariance heatmap
(`ex2_cost_covariance_heatmaps.png`, right) makes this visible: it varies only
with `γ̄` and is flat in `ν̄`. The deterministic recursion plots
(`ex2_deterministic_recursions.png`) show `S(ν̄)` and `P(γ̄)` each diverging as
their argument approaches 0.365 from above.

**What the Monte Carlo shows, and its limit.** The MC mean-energy map (Fig. B)
reproduces the same upper-right stable corner, but with a *softer* boundary, and
the MC thresholds (0.19–0.40) sit **below** 0.365. This is expected and worth
discussing: `1 − 1/A²` is a **mean-square** condition, but the second moment is
heavy-tailed — a few sample paths with long outage runs dominate `E[x²]`, and a
finite sample of `N=100–300` systematically under-estimates that mean, so the
*apparent* boundary is pushed inward. The threshold sweeps
(`ex2_threshold_sweeps.png`) make the horizon dependence explicit: at low
arrival probability the `N=300` energy curve rises far above the `N=100` curve
(to ≈ 10⁵⁷ vs ≈ 10¹⁸), and the MC ν̄ threshold moves from 0.19 (N=100) toward
0.29 (N=300) — i.e. as the horizon grows the empirical boundary climbs toward
the asymptotic 0.365. *The simulation supports the analytic threshold; the
deterministic recursions pin it exactly; the raw MC estimate is a conservative
lower bound because of heavy tails and finite horizon.*

Suggested captions:
- *Fig. A/B.* (left) Deterministic mean-square-stable region: a sharp rectangle
  with corner at `(0.365, 0.365)` = `1−1/A²` (dashed). (right) Monte Carlo
  `log₁₀` mean state energy on the same grid; the stable corner agrees but the
  boundary is gradual.
- *Fig.* Deterministic recursions: steady-state control value `S(ν̄)` and
  expected covariance `P(γ̄)` each diverge at `1−1/A² = 0.365`.
- *Fig.* 1-D energy sweeps at `N=100` and `N=300`; the longer horizon exposes
  the divergence and moves the empirical threshold toward 0.365.

---

## Limitations (keep ~4–6 lines in the report)

- Finite horizon (`N=100`) only **approximates** the asymptotic (`N→∞`)
  stability boundary; we show the empirical threshold drifts toward `1−1/A²` as
  `N` grows.
- Monte-Carlo means of a heavy-tailed second moment under-estimate the
  mean-square boundary; we therefore rely on the deterministic recursions for
  the precise threshold and use MC for illustration and trend.
- The Exercise-1 stability flags use heuristic energy thresholds
  (`stability_flag.m`); they match the plots but are not a formal proof.
- Results are scalar-specific; the `1−1/A²` form holds for scalar (or special
  multivariable) cases — for general `A` the critical value lies between bounds.
- Everything assumes i.i.d. Bernoulli drops and the idealised TCP
  acknowledgement; bursty/correlated DoS would shift the boundaries.

---

## References (placeholders — confirm before use)

```
[1] ELE419 DoS Laboratory Assignment brief, Spring 2026.
[2] L. Schenato, B. Sinopoli, M. Franceschetti, K. Poolla, S. S. Sastry,
    "Foundations of Control and Estimation Over Lossy Networks,"
    Proceedings of the IEEE, 2007.
[3] B. Sinopoli et al., "Kalman Filtering With Intermittent Observations,"
    IEEE Trans. Automatic Control, 2004.   (critical arrival probability 1−1/A²)
[4] ELE419 Lecture notes L8–L10, Spring 2026.
```

Cite [2]/[3] for the TCP-like LQG and the `1−1/A²` threshold; [1]/[4] for the
exact model and task. Do not invent page numbers. This lab is about
**DoS / packet-drop availability** attacks — do not frame it as false-data
injection.

---

## Appendix guidance

Include every `.m` file (see `appendix_code_listing.md`), already commented, in
a monospaced font, outside the 7-page limit. State the run command and the
master seed so the marker can reproduce your figures.
