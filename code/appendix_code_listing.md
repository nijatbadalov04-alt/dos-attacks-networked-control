# Appendix — MATLAB code listing (guidance)

The brief requires **all MATLAB code used for the numerical results, commented**,
in an appendix that does *not* count toward the 7-page limit. Every `.m` file in
this project is already commented (purpose, the *why* behind non-obvious steps,
the indexing convention, and the safe handling of `R=0`). Include them all.

## Recommended order (top-down, easiest to read)

1. `main_ELE419_DoS_lab.m` — driver / configuration
2. `base_params.m` — shared configuration struct
3. `run_exercise1.m` — Exercise 1
4. `run_exercise2.m` — Exercise 2
5. `simulate_tcp_case.m` — single-trajectory simulator
6. `kalman_tcp_update.m` — TCP-like Kalman update
7. `tcp_control_gain.m` — control MARE / steady-state gain
8. `tcp_expected_covariance.m` — expected error-covariance recursion
9. `stability_flag.m` — heuristic Ex-1 stability label
10. `plot_case_results.m`, `make_summary_tables.m`, `save_figure.m` — helpers

State at the top of the appendix: *"Run `main_ELE419_DoS_lab` from a clean
session; master seed = 2026; MATLAB R2025b; no toolboxes required."*

## Build a single combined listing automatically

Keep the appendix in sync with the real files by generating it from them.

**PowerShell** (produces `outputs/code_appendix.txt` in the recommended order):

```powershell
$order = 'main_ELE419_DoS_lab','base_params','run_exercise1','run_exercise2',
         'simulate_tcp_case','kalman_tcp_update','tcp_control_gain',
         'tcp_expected_covariance','stability_flag','plot_case_results',
         'make_summary_tables','save_figure'
$out = 'outputs\code_appendix.txt'
Remove-Item $out -ErrorAction SilentlyContinue
foreach ($f in $order) {
    "%% ===== $f.m =====`r`n" | Out-File $out -Append -Encoding utf8
    Get-Content "$f.m"        | Out-File $out -Append -Encoding utf8
    "`r`n"                    | Out-File $out -Append -Encoding utf8
}
```

**MATLAB** alternative (one nicely formatted HTML/PDF per file via `publish`):

```matlab
files = ["main_ELE419_DoS_lab","base_params","run_exercise1","run_exercise2", ...
         "simulate_tcp_case","kalman_tcp_update","tcp_control_gain", ...
         "tcp_expected_covariance","stability_flag","plot_case_results", ...
         "make_summary_tables","save_figure"];
for f = files
    publish(f + ".m", struct('evalCode',false,'format','pdf', ...
            'outputDir',fullfile('outputs','code_pdf')));
end
```

## Putting it in the report

- **LaTeX:** use the `listings` package and `\lstinputlisting[language=Matlab]{main_ELE419_DoS_lab.m}`
  per file, or `\lstinputlisting{outputs/code_appendix.txt}` for the combined
  file. A small monospaced font (`\footnotesize` / `\small`) keeps it tidy.
- **Word:** paste `code_appendix.txt` into a Consolas/Courier-New, single-spaced
  block; insert a page break before the appendix.

Do not edit the code to "look nicer" for the appendix and then submit different
code that generated the figures — the listing must be exactly what produced
`outputs/`.
