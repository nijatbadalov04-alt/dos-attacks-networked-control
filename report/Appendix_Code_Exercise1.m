%% Appendix A: MATLAB code for Exercise 1
% ELE419 Cybersecurity for Control Systems
% Nijat Badalov (250216213), University of Sheffield, Spring 2026
%
% Denial-of-service packet dropouts on a scalar TCP-like LQG control loop.
% This script reproduces every Exercise 1 result. The six attack scenarios are
% simulated for A = 1.1, N = 100 and process-noise variance 2, using a
% steady-state TCP-like LQG controller together with a TCP-like Kalman filter.
% For each scenario one representative trajectory is plotted, and 300 Monte
% Carlo trials are used for the summary statistics. The parameter values follow
% the ELE419 DoS laboratory brief. Only base MATLAB is used (no toolboxes).
%
% Run from a clean session:  >> Appendix_Code_Exercise1

clear; close all; clc;

masterSeed = 2026;     % fixed seed so the figures and table are reproducible
M          = 300;      % Monte Carlo trials per scenario

figDir = fullfile(pwd, 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% ---- plant and cost (Exercise 1 parameters from the laboratory brief) ------
P.A   = 1.1;   P.B = 1;  P.C = 1;  P.N = 100;   % scalar dynamics and horizon
P.Qxx = 1;     P.Quu = 1;                       % LQ state / input penalties
P.Q   = 2;                                      % process-noise variance Var(w)
P.P0  = 0.01;                                   % initial error covariance
P.x0  = [];                                     % random x0 ~ N(0,P0) each trial
P.validate = true;

% ---- six scenarios: {name, nu_bar, gamma_bar, Var(z)=R, description} -------
Cs = {
    'case1', 1.00, 1.00, 0.0, 'Ideal: no loss, no measurement noise';
    'case2', 1.00, 1.00, 1.0, 'No loss, measurement noise N(0,1)';
    'case3', 0.05, 1.00, 0.5, 'Severe actuation loss, perfect sensing';
    'case4', 1.00, 0.05, 0.5, 'Perfect actuation, severe sensing loss';
    'case5', 0.00, 1.00, 0.1, 'No actuation, sensing available';
    'case6', 1.00, 0.00, 0.1, 'Actuation available, no sensing'};
nC = size(Cs, 1);

rows    = struct([]);
repOuts = cell(1, nC);

for c = 1:nC
    p = P;
    p.nu_bar    = Cs{c,2};
    p.gamma_bar = Cs{c,3};
    p.R         = Cs{c,4};
    [p.L, ~, p.gain_converged] = tcp_control_gain(p);   % steady-state gain once

    % representative trajectory for the per-case figure
    rng(masterSeed + 100*c + 1);
    repOuts{c} = simulate_tcp_case(p);

    % Monte Carlo statistics
    Jr=zeros(1,M); mc=zeros(1,M); rx=zeros(1,M); re=zeros(1,M);
    mx=zeros(1,M); xN2=zeros(1,M); eu=zeros(1,M); eg=zeros(1,M);
    for t = 1:M
        rng(masterSeed + 100*c + t);
        o = simulate_tcp_case(p);
        Jr(t)=o.J_realised; mc(t)=o.mean_cost; rx(t)=o.rms_x; re(t)=o.rms_e;
        mx(t)=o.max_absx;  xN2(t)=o.xN2;   eu(t)=o.emp_nu; eg(t)=o.emp_gamma;
    end
    flag = stability_flag(median(xN2), mean(rx));

    info.name  = Cs{c,1};
    info.title = sprintf('%s: %s   (\\nu=%.2f, \\gamma=%.2f, Var z=%.2g)', ...
                         Cs{c,1}, Cs{c,5}, Cs{c,2}, Cs{c,3}, Cs{c,4});
    plot_case_results(repOuts{c}, info, figDir);

    r.Case=Cs{c,1}; r.nu_bar=Cs{c,2}; r.gamma_bar=Cs{c,3}; r.Var_z=Cs{c,4};
    r.mean_nu=mean(eu); r.mean_gamma=mean(eg);
    r.J_final_mean=mean(Jr); r.mean_cost=mean(mc);
    r.rms_x_mean=mean(rx); r.rms_x_std=std(rx);
    r.rms_e_mean=mean(re); r.rms_e_std=std(re);
    r.max_absx_mean=mean(mx); r.stability=flag; r.note=Cs{c,5};
    if c==1; rows=r; else; rows(c)=r; end

    fprintf(['%s | nu=%.2f gamma=%.2f Var(z)=%.2g | RMSx=%.3g RMSe=%.3g ', ...
             'Jmean=%.3g | %s\n'], Cs{c,1},Cs{c,2},Cs{c,3},Cs{c,4}, ...
             mean(rx),mean(re),mean(Jr),flag);
end

% ---- summary table ---------------------------------------------------------
T = struct2table(rows);
writetable(T, fullfile(figDir, 'exercise1_summary.csv'));
disp(T);

% ---- cross-case comparison figures -----------------------------------------
plot_exercise1_comparison(repOuts, rows, Cs, figDir);

% ---- steady-state vs forward Riccati gain (controller validation) ----------
validate_steady_state_gain(P);     % A = 1.1, nu_bar = 1 representative

fprintf('\nExercise 1 complete. Figures and exercise1_summary.csv saved in: %s\n', figDir);


%% ========================= local functions =============================

function out = simulate_tcp_case(p)
% Simulate one realisation of the scalar TCP-like LQG loop over horizon p.N:
%   x_{k+1} = A x_k + nu_k B u_k + w_k ,  w_k ~ N(0,Q)
%   y_k     = gamma_k C x_k + z_k       ,  z_k ~ N(0,R)
%   u_k     = -L xhat_k     (xhat from the TCP-like Kalman filter)
% Time k = 0..N maps to MATLAB index k+1; inputs/measurements use k = 0..N-1.
% The caller seeds rng before each call, so trials are reproducible.

    if ~isfield(p,'validate') || p.validate
        assert(p.nu_bar>=0 && p.nu_bar<=1,   'nu_bar must be in [0,1]');
        assert(p.gamma_bar>=0 && p.gamma_bar<=1, 'gamma_bar must be in [0,1]');
        assert(p.N>=2 && p.N==round(p.N),    'N must be an integer >= 2');
        assert(p.Q>=0 && p.R>=0,             'noise variances must be >= 0');
    end
    A=p.A; B=p.B; C=p.C; N=p.N; Qxx=p.Qxx; Quu=p.Quu;

    if isfield(p,'L') && ~isempty(p.L)
        L = p.L;
        if isfield(p,'gain_converged'); gc=p.gain_converged; else; gc=NaN; end
    else
        [L,~,gc] = tcp_control_gain(p);
    end

    % all randomness drawn up front (indices 1..N hold the k = 0..N-1 values)
    nu    = double(rand(1,N) < p.nu_bar);     % actuation arrivals
    gamma = double(rand(1,N) < p.gamma_bar);  % sensing arrivals
    w     = sqrt(p.Q)*randn(1,N);             % process noise
    z     = sqrt(p.R)*randn(1,N);             % measurement noise
    if isempty(p.x0); x0 = sqrt(p.P0)*randn; else; x0 = p.x0; end

    x=zeros(1,N+1); x(1)=x0;          % states x_0..x_N
    xhat=zeros(1,N+1);                % posterior estimates xhat_0..xhat_{N-1}
    u=zeros(1,N); y=zeros(1,N);
    Ppost=zeros(1,N+1); Ppost(1)=p.P0;
    Ppred=zeros(1,N+1); Ppred(1)=p.P0;
    Kk=zeros(1,N); mu0=0;

    % k = 0: prior, gated correction with y_0, then u_0
    y(1) = gamma(1)*C*x(1) + z(1);
    s0   = C*p.P0*C' + p.R;
    if gamma(1)==1 && s0>0
        K0=(p.P0*C')/s0; xhat(1)=mu0+K0*(y(1)-C*mu0); Ppost(1)=p.P0-K0*C*p.P0; Kk(1)=K0;
    else
        xhat(1)=mu0; Ppost(1)=p.P0;
    end
    u(1)=-L*xhat(1);
    x(2)=A*x(1)+nu(1)*B*u(1)+w(1);

    % k = 1..N-1: TCP predict/correct, control, propagate
    for k=1:N-1
        j=k+1;
        y(j)=gamma(j)*C*x(j)+z(j);
        [xhat(j),Ppost(j),~,Ppred(j),Kk(j)] = kalman_tcp_update( ...
            p, y(j), xhat(k), u(k), nu(k), gamma(j), Ppost(k));
        u(j)=-L*xhat(j);
        x(j+1)=A*x(j)+nu(j)*B*u(j)+w(j);
    end

    ee = x(1:N)-xhat(1:N);                              % estimation error
    stage = Qxx*x(1:N).^2 + nu.*(Quu*u.^2);             % per-step cost
    J_run = cumsum(stage) + Qxx*x(2:N+1).^2;
    J_run = [Qxx*x(1)^2, J_run];
    J_realised = Qxx*x(N+1)^2 + sum(stage);

    out.p=p; out.k_state=0:N; out.k_io=0:N-1;
    out.x=x; out.xhat=xhat(1:N); out.ee=ee; out.u=u; out.y=y;
    out.nu=nu; out.gamma=gamma; out.Ppost=Ppost(1:N); out.K=Kk;
    out.J_run=J_run; out.L=L; out.gain_converged=gc;
    out.J_realised=J_realised; out.mean_cost=J_realised/N;
    out.emp_nu=mean(nu); out.emp_gamma=mean(gamma);
    out.rms_x=sqrt(mean(x.^2)); out.rms_e=sqrt(mean(ee.^2));
    out.max_absx=max(abs(x)); out.xN2=x(N+1)^2; out.P_final=Ppost(N);
end

function [xhat_upd,P_upd,xhat_pred,P_pred,K] = ...
        kalman_tcp_update(p,y,xhat_prev,u_prev,nu_prev,gamma_cur,P_prev)
% One step of the TCP-like Kalman filter. The previous actuation outcome
% nu_prev is known (TCP acknowledgement), so it enters the prediction
% deterministically; the current sensing arrival gamma_cur gates the
% correction. Division is used instead of inv(); R = 0 is handled safely
% because the innovation covariance C P^- C' + R >= Q > 0 stays positive.
    A=p.A; B=p.B; C=p.C; Q=p.Q; R=p.R;
    xhat_pred = A*xhat_prev + nu_prev*B*u_prev;        % prediction
    P_pred    = A*P_prev*A' + Q;
    innov = C*P_pred*C' + R;
    if gamma_cur==1 && all(innov(:)>0)                 % correction
        K        = (P_pred*C')/innov;
        xhat_upd = xhat_pred + K*(y - C*xhat_pred);
        P_upd    = P_pred - K*C*P_pred;
        P_upd    = 0.5*(P_upd+P_upd.');
    else                                               % no measurement
        K        = zeros(size(P_pred*C'));
        xhat_upd = xhat_pred;
        P_upd    = P_pred;
    end
end

function [L,S,converged] = tcp_control_gain(p)
% Steady-state TCP-like LQG control gain from the modified Riccati recursion
%   L = (B'SB+Quu)^{-1} B'SA ,  S = A'SA + Qxx - nu_bar A'SB L .
% Iterated to its fixed point. For nu_bar below 1 - 1/A^2 the recursion has no
% finite fixed point (S grows without bound) but the gain still converges to
% A/B; converged is then returned false.
    A=p.A; B=p.B; Qxx=p.Qxx; Quu=p.Quu; nu=p.nu_bar;
    S=Qxx; L=0; converged=false; Scap=1e12;
    for it=1:2000
        Ln=(B'*S*B+Quu)\(B'*S*A);
        Sn=A'*S*A+Qxx-nu*(A'*S*B*Ln);
        if abs(Ln-L)<=1e-12*(1+abs(Ln)); L=Ln; S=Sn; converged=true; return; end
        L=Ln; S=Sn;
        if ~isfinite(S)||S>Scap; converged=false; return; end
    end
end

function flag = stability_flag(termEnergy, rmsState)
% Heuristic finite-horizon label from Monte Carlo summaries. A mean-square
% stable loop settles to an O(1-10) energy; an unstable mode inflates the
% terminal energy by many orders of magnitude over the horizon.
    if ~isfinite(termEnergy) || termEnergy>1e3 || rmsState>1e2
        flag='unstable';
    elseif termEnergy>25 || rmsState>10
        flag='marginal';
    else
        flag='stable';
    end
end

function fig = plot_case_results(out, caseInfo, outDir)
% Four-panel diagnostic figure: running cost, state and estimate, estimation
% error with the +/- sqrt(P_k) band, and the packet-arrival sequences.
    N=out.p.N;
    fig=figure('Name',caseInfo.name,'Color','w','Position',[80 60 780 760]);
    subplot(4,1,1);
    plot(out.k_state, log10(max(out.J_run,realmin)),'LineWidth',1.3);
    grid on; xlim([0 N]); xlabel('k'); ylabel('log_{10} J(k)');
    title('Running realised cost');
    subplot(4,1,2);
    stairs(out.k_state,out.x,'LineWidth',1.1); hold on;
    stairs(out.k_io,out.xhat,'--','LineWidth',1.1);
    grid on; xlim([0 N]); xlabel('k'); ylabel('state');
    legend({'x_k','xhat_k'},'Location','best'); title('State and TCP-like estimate');
    subplot(4,1,3);
    band=sqrt(max(out.Ppost,0));
    plot(out.k_io,out.ee,'LineWidth',1.1); hold on;
    plot(out.k_io,band,':','LineWidth',1.0); plot(out.k_io,-band,':','LineWidth',1.0);
    grid on; xlim([0 N]); xlabel('k'); ylabel('x_k - xhat_k');
    title('Estimation error (dotted: \pm one std. dev. \surd P_k)','Interpreter','tex');
    subplot(4,1,4);
    yyaxis left;  stairs(out.k_io,out.gamma,'LineWidth',1.0); ylim([-0.2 1.2]); ylabel('\gamma_k');
    yyaxis right; stairs(out.k_io,out.nu,'LineWidth',1.0);    ylim([-0.2 1.2]); ylabel('\nu_k');
    xlim([0 N]); xlabel('k'); title('Packet arrivals (\gamma_k sensing, \nu_k actuation)');
    if isfield(caseInfo,'title')
        try; sgtitle(caseInfo.title,'FontWeight','bold'); catch; end
    end
    save_fig(fig, outDir, caseInfo.name);
end

function plot_exercise1_comparison(repOuts, results, Cs, outDir)
% Cross-case overlays: running cost and RMS state by scenario.
    nC=numel(repOuts); names=Cs(:,1);
    f1=figure('Color','w','Position',[80 80 760 440]); cols=lines(nC); hold on;
    for c=1:nC
        plot(repOuts{c}.k_state, log10(max(repOuts{c}.J_run,realmin)), ...
             'LineWidth',1.4,'Color',cols(c,:));
    end
    grid on; xlabel('k'); ylabel('log_{10} J(k)'); legend(names,'Location','northwest');
    title('Exercise 1: running realised cost (representative run per case)');
    save_fig(f1,outDir,'ex1_cost_overlay');

    f2=figure('Color','w','Position',[80 80 760 440]);
    rx=[results.rms_x_mean]; sx=[results.rms_x_std];
    bar(rx); hold on; lowErr=min(sx,rx*0.999);
    errorbar(1:nC,rx,lowErr,sx,'k','LineStyle','none','CapSize',8);
    set(gca,'XTick',1:nC,'XTickLabel',names,'YScale','log');
    ylabel('RMS state (MC mean \pm std)'); grid on;
    title('Exercise 1: RMS state by case');
    save_fig(f2,outDir,'ex1_rms_state_bar');
end

function validate_steady_state_gain(P)
% Compare the steady-state control gain with the forward Riccati-style
% recursion of the provided implementation (nu_bar = 1, A = 1.1). Reports the
% step at which the forward gain matches the steady-state value and the
% residual deviation thereafter.
    p=P; p.nu_bar=1;
    [Lss,~,~]=tcp_control_gain(p);
    A=p.A; B=p.B; Qxx=p.Qxx; Quu=p.Quu; N=p.N;
    S=Qxx; Lk=zeros(1,N);
    for k=2:N
        Lk(k)=(B'*S*B+Quu)\(B'*S*A);
        S=A'*S*A+Qxx-p.nu_bar*(A'*S*B*Lk(k));
    end
    dev=abs(Lk(2:N)-Lss);
    kConv=find(dev<1e-6,1,'first')+1;       % first step within 1e-6 of L_ss
    fprintf(['\nController check (A=%.3f, nu_bar=1): steady-state gain L=%.4f; ', ...
             'forward Riccati gain matches it to 1e-6 by step %d; max deviation ', ...
             'after step 5 = %.2e.\n'], A, Lss, kConv, max(dev(5:end)));
end

function save_fig(figHandle, outDir, baseName)
% Save a figure as PNG using a light (print-friendly) theme.
    if ~exist(outDir,'dir'); mkdir(outDir); end
    try; theme(figHandle,'light'); catch; end
    set(figHandle,'Color','w');
    ax=findall(figHandle,'Type','axes');
    for a=reshape(ax,1,[]); try; a.Toolbar.Visible='off'; catch; end; end
    png=fullfile(outDir,[baseName '.png']);
    try; exportgraphics(figHandle,png,'Resolution',150);
    catch; try; print(figHandle,png,'-dpng','-r150'); catch; saveas(figHandle,png); end; end
end
