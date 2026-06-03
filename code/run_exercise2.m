function ex2 = run_exercise2(cfg, outDirs)
%RUN_EXERCISE2  Stability search over (nu_bar, gamma_bar) for the unstable
%   A = 1.255 scalar system (brief Exercise 2).
%
%   System: A=1.255, B=1, C=1, N=100, w_k ~ N(0,0.5), z_k ~ N(0,0.15).
%   Goal  : find the (nu_bar, gamma_bar) values for which the closed loop is
%           (mean-square) stable, and compare with the analytic critical
%           arrival probability 1 - 1/A^2 = 0.3651.
%
%   TWO COMPLEMENTARY ANALYSES.
%   (A) DETERMINISTIC (precise). The mean-square stability boundary is exactly
%       where the two coupled recursions lose boundedness:
%         * control:   the MARE  S <- A'SA + Qxx - nu_bar A'SB(B'SB+Quu)^{-1}B'SA
%                      has a finite fixed point iff nu_bar > 1 - 1/A^2;
%         * estimator: the expected-covariance recursion (tcp_expected_covariance)
%                      is bounded iff gamma_bar > 1 - 1/A^2.
%       In the TCP-like setup the estimator recursion does NOT depend on nu and
%       the control recursion does NOT depend on gamma (separation), so the
%       stable region is the rectangle {nu_bar > 0.365} x {gamma_bar > 0.365}.
%       We locate both thresholds by bisection on recursion boundedness.
%   (B) MONTE CARLO (illustrative). We simulate the full stochastic loop on a
%       grid and a fine 1-D sweep and report the mean state energy. This shows
%       the trend and confirms the boundary, but a finite-sample / finite-horizon
%       mean systematically UNDER-estimates the mean-square boundary because the
%       second moment is heavy-tailed; we make that limitation explicit.
%
%   cfg     : struct with .masterSeed, .ex2grid, .ex2trials, .ex2horizon.
%   outDirs : struct with fields ex2, figures, tables.

    fprintf('\n==== EXERCISE 2  (A=1.255, N=%d, Var(w)=0.5, Var(z)=0.15) ====\n', cfg.ex2horizon);

    p0          = base_params();
    p0.A        = 1.255;
    p0.Q        = 0.5;                 % process-noise variance
    p0.R        = 0.15;                % measurement-noise variance
    p0.N        = cfg.ex2horizon;      % validation horizon (brief value 100)
    p0.validate = false;               % skip per-call asserts in the sweeps
    crit        = 1 - 1/p0.A^2;        % analytic critical arrival probability

    g  = cfg.ex2grid(:).';
    nb = g; gb = g;
    nN = numel(nb); nG = numel(gb);

    % ====================================================================
    % (A) DETERMINISTIC mean-square stability boundary
    % ====================================================================
    nu_crit_det  = bisect_threshold(@(v) ctrl_bounded(setfield(p0,'nu_bar',v)),    0, 1, 1e-3); %#ok<SFLD>
    gam_crit_det = bisect_threshold(@(v) cov_bounded( setfield(p0,'gamma_bar',v)), 0, 1, 1e-3); %#ok<SFLD>

    % deterministic stable map on the grid (outer product of the two channels)
    ctrlOK = false(1,nN); for i = 1:nN; ctrlOK(i) = ctrl_bounded(setfield(p0,'nu_bar',nb(i))); end %#ok<SFLD>
    covOK  = false(1,nG); for j = 1:nG; covOK(j)  = cov_bounded( setfield(p0,'gamma_bar',gb(j))); end %#ok<SFLD>
    stableDet = double(covOK(:) * ctrlOK(:).');     % nG x nN, 1 = MS-stable

    % deterministic 1-D curves: steady control value S(nu) and covariance P(gamma)
    gf   = 0:0.01:1;
    Sdet = nan(size(gf)); Pdet = nan(size(gf));
    for i = 1:numel(gf)
        [Si, ci] = ctrl_steady(setfield(p0,'nu_bar',gf(i)));            %#ok<SFLD>
        if ci; Sdet(i) = Si; else; Sdet(i) = NaN; end                  % NaN => diverges
        [Pi, bi] = tcp_expected_covariance(setfield(p0,'gamma_bar',gf(i))); %#ok<SFLD>
        if bi; Pdet(i) = Pi; else; Pdet(i) = NaN; end
    end

    % ====================================================================
    % (B) MONTE CARLO grid: mean state energy, cost, final covariance
    % ====================================================================
    M = cfg.ex2trials;
    logEnergy = nan(nG, nN);
    meanCost  = nan(nG, nN);
    finalCov  = nan(nG, nN);

    t0 = tic;
    for in = 1:nN
        p = p0; p.nu_bar = nb(in);
        [p.L, ~, p.gain_converged] = tcp_control_gain(p);   % gain depends on nu only
        for ig = 1:nG
            p.gamma_bar = gb(ig);
            sumX2 = zeros(1, p.N+1); eJ = zeros(1,M); eP = zeros(1,M);
            for t = 1:M
                rng(cfg.masterSeed + 5000 + 137*in + t);
                o     = simulate_tcp_case(p);
                sumX2 = sumX2 + o.x.^2;
                eJ(t) = o.mean_cost;
                eP(t) = o.P_final;
            end
            V = sumX2 / M;
            logEnergy(ig,in) = log10(max(mean(V), realmin));
            meanCost(ig,in)  = mean(eJ);
            finalCov(ig,in)  = mean(eP);
        end
    end
    fprintf('  MC grid: %d x %d points, %d trials each (%.1f s).\n', nN, nG, M, toc(t0));

    % MC 1-D sweeps (fine grid, two horizons) --------------------------------
    Ns        = [cfg.ex2horizon, max(3*cfg.ex2horizon, 300)];
    sweep_nu  = run_sweep(p0, gf, 'nu',    1.0, Ns, cfg);
    sweep_gam = run_sweep(p0, gf, 'gamma', 1.0, Ns, cfg);

    % ====================================================================
    % FIGURES
    % ====================================================================
    % (1) deterministic stable region  +  MC mean state energy
    f1 = figure('Color','w','Position',[60 60 1000 430]);
    subplot(1,2,1);
    draw_heatmap(nb, gb, stableDet, crit, '', false);
    title('(A) Deterministic MS-stable region');
    subplot(1,2,2);
    draw_heatmap(nb, gb, logEnergy, crit, 'log_{10} mean state energy', true);
    title('(B) Monte Carlo log_{10} E[(1/N)\Sigma x_k^2]');
    save_figure(f1, outDirs.ex2, 'ex2_stability_heatmaps');

    % (2) MC cost map  +  MC final covariance (depends on gamma only)
    f2 = figure('Color','w','Position',[60 60 1000 430]);
    subplot(1,2,1);
    draw_heatmap(nb, gb, log10(max(meanCost,realmin)), crit, 'log_{10} mean cost', true);
    title('Cost map: log_{10} mean per-step cost');
    subplot(1,2,2);
    draw_heatmap(nb, gb, log10(max(finalCov,realmin)), crit, 'log_{10} final covariance', true);
    title('Estimator: log_{10} final P (varies with \gamma only)');
    save_figure(f2, outDirs.ex2, 'ex2_cost_covariance_heatmaps');

    % (3) MC threshold sweeps (energy vs arrival prob, two horizons)
    f3 = figure('Color','w','Position',[60 60 1000 430]);
    subplot(1,2,1);
    plot_sweep(sweep_nu, crit, '\nu_{bar} (\gamma_{bar}=1)');
    title('Vary actuation arrival (MC)');
    subplot(1,2,2);
    plot_sweep(sweep_gam, crit, '\gamma_{bar} (\nu_{bar}=1)');
    title('Vary sensing arrival (MC)');
    save_figure(f3, outDirs.ex2, 'ex2_threshold_sweeps');

    % (4) deterministic recursions: S(nu) and Pbar(gamma) diverge at 1-1/A^2
    f4 = figure('Color','w','Position',[60 60 1000 430]);
    subplot(1,2,1);
    semilogy(gf, Sdet, 'LineWidth', 1.6); hold on; yl = ylim;
    plot([crit crit], yl, 'k--', 'LineWidth', 1.2);
    grid on; xlabel('\nu_{bar}'); ylabel('steady-state control value S');
    legend({'S(\nu_{bar})', sprintf('1-1/A^2=%.3f',crit)}, 'Location','northeast');
    title('Control MARE diverges as \nu_{bar}\downarrow 1-1/A^2');
    subplot(1,2,2);
    semilogy(gf, Pdet, 'LineWidth', 1.6); hold on; yl = ylim;
    plot([crit crit], yl, 'k--', 'LineWidth', 1.2);
    grid on; xlabel('\gamma_{bar}'); ylabel('steady-state expected covariance P');
    legend({'P(\gamma_{bar})', sprintf('1-1/A^2=%.3f',crit)}, 'Location','northeast');
    title('Expected covariance diverges as \gamma_{bar}\downarrow 1-1/A^2');
    save_figure(f4, outDirs.ex2, 'ex2_deterministic_recursions');

    % ====================================================================
    % SUMMARY (console + table)
    % ====================================================================
    fprintf('\n  --- Exercise 2 stability summary ---\n');
    fprintf('  Analytic critical arrival probability   1 - 1/A^2     = %.4f\n', crit);
    fprintf('  Deterministic nu threshold  (control MARE bounded)    = %.4f\n', nu_crit_det);
    fprintf('  Deterministic gamma threshold (expected cov. bounded) = %.4f\n', gam_crit_det);
    fprintf('  MC fine-sweep nu threshold   : short N=%d ~ %.3f , long N=%d ~ %.3f\n', ...
            Ns(1), sweep_nu.thr_short, Ns(2), sweep_nu.thr_long);
    fprintf('  MC fine-sweep gamma threshold: short N=%d ~ %.3f , long N=%d ~ %.3f\n', ...
            Ns(1), sweep_gam.thr_short, Ns(2), sweep_gam.thr_long);
    fprintf('  => STABLE iff  nu_bar > %.3f  AND  gamma_bar > %.3f  (TCP separation).\n', ...
            nu_crit_det, gam_crit_det);
    fprintf('     (MC means under-estimate the boundary: heavy-tailed 2nd moment + finite N.)\n');

    rows = struct();
    rows(1).quantity = 'analytic_1_minus_1overA2';     rows(1).value = crit;
    rows(2).quantity = 'det_nu_threshold_MARE';        rows(2).value = nu_crit_det;
    rows(3).quantity = 'det_gamma_threshold_cov';      rows(3).value = gam_crit_det;
    rows(4).quantity = 'mc_nu_threshold_shortN';       rows(4).value = sweep_nu.thr_short;
    rows(5).quantity = 'mc_nu_threshold_longN';        rows(5).value = sweep_nu.thr_long;
    rows(6).quantity = 'mc_gamma_threshold_shortN';    rows(6).value = sweep_gam.thr_short;
    rows(7).quantity = 'mc_gamma_threshold_longN';     rows(7).value = sweep_gam.thr_long;
    make_summary_tables(rows, outDirs.tables, 'exercise2_thresholds');

    ex2 = struct('nb',nb,'gb',gb,'logEnergy',logEnergy,'meanCost',meanCost, ...
                 'finalCov',finalCov,'stableDet',stableDet,'crit',crit, ...
                 'nu_crit_det',nu_crit_det,'gam_crit_det',gam_crit_det, ...
                 'gf',gf,'Sdet',Sdet,'Pdet',Pdet, ...
                 'sweep_nu',sweep_nu,'sweep_gam',sweep_gam);
    save(fullfile(outDirs.ex2,'exercise2_grid.mat'),'-struct','ex2');
end

% ========================================================================
%  DETERMINISTIC boundedness tests (mean-square stability of the recursions)
% ========================================================================
function [Sss, bounded] = ctrl_steady(p)
%CTRL_STEADY  Steady-state control value S from the control MARE, iterated long
%   enough to resolve the stability threshold (the gain helper tcp_control_gain
%   uses a short cap suited to the gain, not to threshold detection). Returns
%   the converged S and a boundedness flag (false => no finite fixed point).
    A = p.A; B = p.B; Qxx = p.Qxx; Quu = p.Quu; nu = p.nu_bar;
    S = Qxx; cap = 1e8; bounded = true;
    for it = 1:1e5
        L  = (B'*S*B + Quu) \ (B'*S*A);
        Sn = A'*S*A + Qxx - nu*(A'*S*B*L);
        if abs(Sn - S) <= 1e-12*(1 + abs(Sn)); Sss = Sn; bounded = true; return; end
        S = Sn;
        if ~isfinite(S) || S > cap; Sss = S; bounded = false; return; end
    end
    Sss = S; bounded = isfinite(S) && S <= cap;
end

function tf = ctrl_bounded(p)
%CTRL_BOUNDED  True iff the control MARE has a finite fixed point, i.e. a
%   stabilising steady-state controller exists (nu_bar above the critical
%   arrival probability).
    [~, tf] = ctrl_steady(p);
end

function tf = cov_bounded(p)
%COV_BOUNDED  True iff the expected error-covariance recursion stays bounded.
    [~, tf] = tcp_expected_covariance(p, 1e5, 1e8);
end

function thr = bisect_threshold(testFun, lo, hi, tol)
%BISECT_THRESHOLD  Smallest value in [lo,hi] for which the (monotone) testFun
%   returns true, located by bisection to resolution tol.
    if ~testFun(hi); thr = NaN; return; end     % never stable
    if  testFun(lo); thr = lo;  return; end      % always stable
    while hi - lo > tol
        mid = 0.5*(lo + hi);
        if testFun(mid); hi = mid; else; lo = mid; end
    end
    thr = hi;
end

% ========================================================================
%  MONTE CARLO helpers
% ========================================================================
function s = run_sweep(p0, gf, which, fixedVal, Ns, cfg)
%RUN_SWEEP  1-D Monte Carlo sweep of one arrival probability at two horizons.
%   Returns log mean state energy and a (finite-horizon) MS-growth threshold.
    M    = max(150, cfg.ex2trials);
    logE = nan(numel(Ns), numel(gf));
    grow = nan(numel(Ns), numel(gf));
    for hi = 1:numel(Ns)
        p = p0; p.N = Ns(hi);
        for i = 1:numel(gf)
            if strcmp(which,'nu')
                p.nu_bar = gf(i); p.gamma_bar = fixedVal;
            else
                p.gamma_bar = gf(i); p.nu_bar = fixedVal;
            end
            [p.L,~,p.gain_converged] = tcp_control_gain(p);
            sumX2 = zeros(1, p.N+1);
            for t = 1:M
                rng(cfg.masterSeed + 9000 + 311*hi + 17*i + t);
                o = simulate_tcp_case(p);
                sumX2 = sumX2 + o.x.^2;
            end
            V = sumX2 / M;
            logE(hi,i) = log10(max(mean(V), realmin));
            midI = round(p.N/2)+1;
            grow(hi,i) = V(end) / max(V(midI), realmin);   % MS growth ratio
        end
    end
    s.gf        = gf;
    s.Ns        = Ns;
    s.logE      = logE;
    s.grow      = grow;
    s.thr_short = first_stable(gf, grow(1,:) < 3);
    s.thr_long  = first_stable(gf, grow(2,:) < 3);
    s.which     = which;
end

function thr = first_stable(grid, stableRow)
%FIRST_STABLE  Smallest grid value at which the row becomes (and stays) stable.
    thr = NaN;
    for i = 1:numel(grid)
        if stableRow(i) && all(stableRow(i:end)); thr = grid(i); return; end
    end
end

function draw_heatmap(nb, gb, Z, crit, cbarLabel, useColorbar)
%DRAW_HEATMAP  imagesc heatmap with analytic critical-threshold guide lines.
    imagesc(nb, gb, Z); set(gca,'YDir','normal'); hold on;
    xlabel('\nu_{bar} (actuation arrival)');
    ylabel('\gamma_{bar} (sensing arrival)');
    if useColorbar; cb = colorbar; cb.Label.String = cbarLabel; end
    plot([crit crit], [gb(1) gb(end)], 'w--', 'LineWidth', 1.5);   % nu-critical
    plot([nb(1) nb(end)], [crit crit], 'w--', 'LineWidth', 1.5);   % gamma-critical
    axis([nb(1) nb(end) gb(1) gb(end)]);
end

function plot_sweep(s, crit, xlab)
%PLOT_SWEEP  log mean state energy vs arrival probability for both horizons.
    plot(s.gf, s.logE(1,:), 'LineWidth', 1.5); hold on;
    plot(s.gf, s.logE(2,:), 'LineWidth', 1.5);
    yl = ylim;
    plot([crit crit], yl, 'k--', 'LineWidth', 1.2);
    grid on; xlabel(xlab); ylabel('log_{10} mean state energy');
    legend({sprintf('N=%d',s.Ns(1)), sprintf('N=%d',s.Ns(2)), ...
            sprintf('1-1/A^2=%.3f',crit)}, 'Location','northeast');
end
