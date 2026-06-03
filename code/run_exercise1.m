function results = run_exercise1(cfg, outDirs)
%RUN_EXERCISE1  Six DoS cases for the A = 1.1 scalar system (brief Exercise 1).
%
%   System: A=1.1, B=1, C=1, N=100, w_k ~ N(0,2). Each of the six cases sets a
%   different (nu_bar, gamma_bar, Var(z_k)) per the brief. For every case we
%       * compute the steady-state TCP-like control gain once,
%       * draw one reproducible REPRESENTATIVE trajectory for the plots,
%       * run cfg.mcTrials Monte Carlo trials for mean/std metrics,
%       * save a four-panel figure and a summary-table row.
%
%   cfg     : struct with .masterSeed and .mcTrials.
%   outDirs : struct with fields ex1, figures, tables.
%   results : struct array of per-case summary rows (also written to disk).

    fprintf('\n==== EXERCISE 1  (A=1.1, N=100, Var(w)=2) ====\n');

    p0     = base_params();
    p0.A   = 1.1;
    p0.N   = 100;
    p0.Q   = 2;                       % process-noise variance Var(w_k) = 2
    nu_crit = 1 - 1/p0.A^2;           % theoretical critical arrival prob (A=1.1)

    % case table: { name, nu_bar, gamma_bar, Var(z)=R, description } -----------
    Cs = {
        'case1', 1.00, 1.00, 0.0, 'Ideal: no loss, no measurement noise';
        'case2', 1.00, 1.00, 1.0, 'No loss, measurement noise N(0,1)';
        'case3', 0.05, 1.00, 0.5, 'Severe actuation loss, perfect sensing';
        'case4', 1.00, 0.05, 0.5, 'Perfect actuation, severe sensing loss';
        'case5', 0.00, 1.00, 0.1, 'No actuation, sensing available';
        'case6', 1.00, 0.00, 0.1, 'Actuation available, no sensing';
        };
    nC = size(Cs, 1);

    repOuts = cell(1, nC);
    for c = 1:nC
        p           = p0;
        p.nu_bar    = Cs{c,2};
        p.gamma_bar = Cs{c,3};
        p.R         = Cs{c,4};                       % measurement-noise variance
        [p.L, ~, p.gain_converged] = tcp_control_gain(p);   % gain once per case

        % ---- representative trajectory (fixed seed -> reproducible plot) ----
        rng(cfg.masterSeed + 100*c + 1);
        repOut       = simulate_tcp_case(p);
        repOuts{c}   = repOut;

        % ---- Monte Carlo over cfg.mcTrials trials ---------------------------
        M  = cfg.mcTrials;
        Jr = zeros(1,M); mc = zeros(1,M); rx = zeros(1,M); re = zeros(1,M);
        mx = zeros(1,M); xN2 = zeros(1,M); eg = zeros(1,M); eu = zeros(1,M);
        for t = 1:M
            rng(cfg.masterSeed + 100*c + t);          % per-trial reproducible seed
            o      = simulate_tcp_case(p);
            Jr(t)  = o.J_realised;  mc(t)  = o.mean_cost; rx(t) = o.rms_x;
            re(t)  = o.rms_e;       mx(t)  = o.max_absx;  xN2(t) = o.xN2;
            eg(t)  = o.emp_gamma;   eu(t)  = o.emp_nu;
        end
        flag = stability_flag(median(xN2), mean(rx));

        % ---- per-case diagnostic figure -------------------------------------
        info.name  = Cs{c,1};
        info.title = sprintf('%s: %s   (\\nu=%.2f, \\gamma=%.2f, Var z=%.2g)', ...
                             Cs{c,1}, Cs{c,5}, Cs{c,2}, Cs{c,3}, Cs{c,4});
        plot_case_results(repOut, info, outDirs.ex1);

        % ---- summary row (MC mean +/- std) ----------------------------------
        r = struct();
        r.Case          = Cs{c,1};
        r.nu_bar        = Cs{c,2};
        r.gamma_bar     = Cs{c,3};
        r.Var_z         = Cs{c,4};
        r.mean_nu       = mean(eu);          % realised mean actuation arrival
        r.mean_gamma    = mean(eg);          % realised mean sensing arrival
        r.J_final_mean  = mean(Jr);
        r.J_final_std   = std(Jr);
        r.mean_cost     = mean(mc);
        r.rms_x_mean    = mean(rx);
        r.rms_x_std     = std(rx);
        r.rms_e_mean    = mean(re);
        r.rms_e_std     = std(re);
        r.max_absx_mean = mean(mx);
        r.gain_conv     = double(repOut.gain_converged);
        r.stability     = flag;
        r.note          = Cs{c,5};

        if c == 1; results = r; else; results(c) = r; end %#ok<AGROW>

        fprintf(['  %s | nu=%.2f gamma=%.2f Var(z)=%.2g | ', ...
                 'RMSx=%.3g RMSe=%.3g Jfinal=%.3g | %s\n'], ...
                 Cs{c,1}, Cs{c,2}, Cs{c,3}, Cs{c,4}, ...
                 mean(rx), mean(re), mean(Jr), flag);
    end

    % ---- export summary table and comparison figures -----------------------
    make_summary_tables(results, outDirs.tables, 'exercise1_summary');
    plot_exercise1_comparison(repOuts, results, Cs, outDirs.figures);

    fprintf('  Theoretical critical arrival probability for A=1.1: 1-1/A^2 = %.3f\n', nu_crit);
    fprintf('  (Cases 3-6 all push one channel below this threshold.)\n');
end

% ========================================================================
function plot_exercise1_comparison(repOuts, results, Cs, outDir)
%PLOT_EXERCISE1_COMPARISON  Cross-case overlays: running cost and RMS state.
    nC   = numel(repOuts);
    names = Cs(:,1);

    % (a) running cost overlay (representative runs) -------------------------
    f1   = figure('Color','w','Position',[80 80 760 440]);
    cols = lines(nC); hold on;
    for c = 1:nC
        plot(repOuts{c}.k_state, log10(max(repOuts{c}.J_run, realmin)), ...
             'LineWidth', 1.4, 'Color', cols(c,:));
    end
    grid on; xlabel('k'); ylabel('log_{10} J(k)');
    legend(names, 'Location', 'northwest');
    title('Exercise 1: running realised cost (representative run per case)');
    save_figure(f1, outDir, 'ex1_cost_overlay');

    % (b) RMS state bar chart with MC std (log y-axis) -----------------------
    f2 = figure('Color','w','Position',[80 80 760 440]);
    rx = [results.rms_x_mean];
    sx = [results.rms_x_std];
    bar(rx); hold on;
    lowErr = min(sx, rx*0.999);                 % keep error bars positive on log axis
    errorbar(1:nC, rx, lowErr, sx, 'k', 'LineStyle', 'none', 'CapSize', 8);
    set(gca, 'XTick', 1:nC, 'XTickLabel', names, 'YScale', 'log');
    ylabel('RMS state (MC mean \pm std)'); grid on;
    title('Exercise 1: RMS state by case');
    save_figure(f2, outDir, 'ex1_rms_state_bar');
end
