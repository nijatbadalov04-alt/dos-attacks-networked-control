%% Appendix B: MATLAB code for Exercise 2
% ELE419 Cybersecurity for Control Systems
% Nijat Badalov (250216213), University of Sheffield, Spring 2026
%
% Mean-square stability boundary of the scalar TCP-like LQG loop under
% denial-of-service packet dropouts, for A = 1.255, N = 100, process-noise
% variance 0.5 and measurement-noise variance 0.15 (ELE419 DoS laboratory
% brief). The stability region in (nu_bar, gamma_bar) is located two ways:
%   (i)  deterministically, from the boundedness of the control Riccati
%        recursion (in nu_bar) and the expected error-covariance recursion
%        (in gamma_bar); these give the exact mean-square boundary;
%   (ii) by Monte Carlo simulation of the full stochastic loop on a grid and
%        on fine one-dimensional sweeps, which supports the trend but is less
%        sharp near the threshold. Only base MATLAB is used (no toolboxes).
%
% Run from a clean session:  >> Appendix_Code_Exercise2

clear; close all; clc;

masterSeed = 2026;
figDir = fullfile(pwd, 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% ---- plant and cost (Exercise 2 parameters) --------------------------------
P.A=1.255; P.B=1; P.C=1; P.N=100; P.Qxx=1; P.Quu=1;
P.Q=0.5;        % process-noise variance Var(w)
P.R=0.15;       % measurement-noise variance Var(z)
P.P0=0.01; P.x0=[]; P.validate=false;
crit = 1 - 1/P.A^2;                 % analytic critical arrival probability

ex2grid    = 0:0.05:1;              % (nu_bar, gamma_bar) grid for the heatmaps
ex2trials  = 80;                    % Monte Carlo trials per grid point
nb = ex2grid; gb = ex2grid; nN=numel(nb); nG=numel(gb);

% ===================================================================
% (i) DETERMINISTIC mean-square stability boundary
% ===================================================================
nu_crit_det  = bisect_threshold(@(v) ctrl_bounded(setfield(P,'nu_bar',v)),    0,1,1e-3); %#ok<SFLD>
gam_crit_det = bisect_threshold(@(v) cov_bounded( setfield(P,'gamma_bar',v)), 0,1,1e-3); %#ok<SFLD>

ctrlOK=false(1,nN); for i=1:nN; ctrlOK(i)=ctrl_bounded(setfield(P,'nu_bar',nb(i))); end %#ok<SFLD>
covOK =false(1,nG); for j=1:nG; covOK(j) =cov_bounded( setfield(P,'gamma_bar',gb(j))); end %#ok<SFLD>
stableDet = double(covOK(:)*ctrlOK(:).');           % nG x nN, 1 = MS-stable

gf=0:0.01:1; Sdet=nan(size(gf)); Pdet=nan(size(gf));
for i=1:numel(gf)
    [Si,ci]=ctrl_steady(setfield(P,'nu_bar',gf(i)));            %#ok<SFLD>
    if ci; Sdet(i)=Si; end
    [Pi,bi]=tcp_expected_covariance(setfield(P,'gamma_bar',gf(i))); %#ok<SFLD>
    if bi; Pdet(i)=Pi; end
end

% ===================================================================
% (ii) MONTE CARLO grid: mean state energy, cost, final covariance
% ===================================================================
M=ex2trials; logEnergy=nan(nG,nN); meanCost=nan(nG,nN); finalCov=nan(nG,nN);
t0=tic;
for in=1:nN
    p=P; p.nu_bar=nb(in);
    [p.L,~,p.gain_converged]=tcp_control_gain(p);   % gain depends on nu_bar only
    for ig=1:nG
        p.gamma_bar=gb(ig);
        sumX2=zeros(1,p.N+1); eJ=zeros(1,M); eP=zeros(1,M);
        for t=1:M
            rng(masterSeed+5000+137*in+t);
            o=simulate_tcp_case(p);
            sumX2=sumX2+o.x.^2; eJ(t)=o.mean_cost; eP(t)=o.P_final;
        end
        V=sumX2/M;
        logEnergy(ig,in)=log10(max(mean(V),realmin));
        meanCost(ig,in)=mean(eJ); finalCov(ig,in)=mean(eP);
    end
end
fprintf('MC grid: %d x %d points, %d trials each (%.1f s).\n', nN,nG,M,toc(t0));

% fine 1-D sweeps at two horizons
Ns=[P.N, 300];
sweep_nu  = run_sweep(P, gf, 'nu',    1.0, Ns, masterSeed);
sweep_gam = run_sweep(P, gf, 'gamma', 1.0, Ns, masterSeed);

% ===================================================================
% FIGURES
% ===================================================================
f1=figure('Color','w','Position',[60 60 1000 430]);
subplot(1,2,1); draw_heatmap(nb,gb,stableDet,crit,'',false);
title('(A) Deterministic MS-stable region');
subplot(1,2,2); draw_heatmap(nb,gb,logEnergy,crit,'log_{10} mean state energy',true);
title('(B) Monte Carlo log_{10} E[(1/N)\Sigma x_k^2]');
save_fig(f1,figDir,'ex2_stability_heatmaps');

f2=figure('Color','w','Position',[60 60 1000 430]);
subplot(1,2,1); draw_heatmap(nb,gb,log10(max(meanCost,realmin)),crit,'log_{10} mean cost',true);
title('Cost map: log_{10} mean per-step cost');
subplot(1,2,2); draw_heatmap(nb,gb,log10(max(finalCov,realmin)),crit,'log_{10} final covariance',true);
title('Estimator: log_{10} final P (varies with \gamma only)');
save_fig(f2,figDir,'ex2_cost_covariance_heatmaps');

f3=figure('Color','w','Position',[60 60 1000 430]);
subplot(1,2,1); plot_sweep(sweep_nu,crit,'\nu_{bar} (\gamma_{bar}=1)'); title('Vary actuation arrival (MC)');
subplot(1,2,2); plot_sweep(sweep_gam,crit,'\gamma_{bar} (\nu_{bar}=1)'); title('Vary sensing arrival (MC)');
save_fig(f3,figDir,'ex2_threshold_sweeps');

f4=figure('Color','w','Position',[60 60 1000 430]);
subplot(1,2,1); semilogy(gf,Sdet,'LineWidth',1.6); hold on; yl=ylim;
plot([crit crit],yl,'k--','LineWidth',1.2); grid on;
xlabel('\nu_{bar}'); ylabel('steady-state control value S');
legend({'S(\nu_{bar})',sprintf('1-1/A^2=%.3f',crit)},'Location','northeast');
title('Control recursion diverges as \nu_{bar}\downarrow 1-1/A^2');
subplot(1,2,2); semilogy(gf,Pdet,'LineWidth',1.6); hold on; yl=ylim;
plot([crit crit],yl,'k--','LineWidth',1.2); grid on;
xlabel('\gamma_{bar}'); ylabel('steady-state expected covariance P');
legend({'P(\gamma_{bar})',sprintf('1-1/A^2=%.3f',crit)},'Location','northeast');
title('Expected covariance diverges as \gamma_{bar}\downarrow 1-1/A^2');
save_fig(f4,figDir,'ex2_deterministic_recursions');

% supporting 3D stability landscape: log10 mean state energy as a surface, with
% a base contour (surfc) and thin dashed threshold lines on the base plane.
% The surface is lightly smoothed and interpolated for DISPLAY ONLY; this does
% not affect any reported threshold, which comes from the recursions above.
[NBb,GBb]=meshgrid(nb,gb);
Zs=smooth2(logEnergy);                                   % light 3x3 mean
nbf=linspace(0,1,101); gbf=linspace(0,1,101); [NBf,GBf]=meshgrid(nbf,gbf);
Zf=interp2(NBb,GBb,Zs,NBf,GBf,'spline');                 % smooth display surface
f5=figure('Color','w','Position',[80 80 780 560]);
surfc(nbf,gbf,Zf,'EdgeColor','none'); shading interp; hold on;
colormap(gca,seqmap); clim([pct(Zf,3) pct(Zf,97)]);      % percentile colour clip
cb=colorbar; cb.Label.String='log_{10} mean state energy';
zb=min(Zf(:))-0.5; zlim([zb max(Zf(:))]);
plot3([crit crit],[0 1],[zb zb],'k--','LineWidth',1.4);  % nu_bar = 0.365 (base)
plot3([0 1],[crit crit],[zb zb],'k--','LineWidth',1.4);  % gamma_bar = 0.365 (base)
xlabel('\nu_{bar}  (actuation arrival prob.)');
ylabel('\gamma_{bar}  (sensing arrival prob.)');
zlabel('log_{10} mean state energy');
title(sprintf('3D stability landscape (dashed thresholds at \\nu_{bar},\\gamma_{bar}\\approx%.3f)',crit));
view(-38,26); grid on; box on; pbaspect([1 1 0.6]);
save_fig(f5,figDir,'ex2_stability_3d');

% ===================================================================
% SUMMARY (console + CSV)
% ===================================================================
fprintf('\nAnalytic critical arrival probability 1-1/A^2 = %.4f\n', crit);
fprintf('Deterministic nu threshold  (control recursion)      = %.4f\n', nu_crit_det);
fprintf('Deterministic gamma threshold (covariance recursion) = %.4f\n', gam_crit_det);
fprintf('Monte Carlo nu threshold    : N=%d ~ %.3f , N=%d ~ %.3f\n', ...
        Ns(1),sweep_nu.thr_short,Ns(2),sweep_nu.thr_long);
fprintf('Monte Carlo gamma threshold : N=%d ~ %.3f , N=%d ~ %.3f\n', ...
        Ns(1),sweep_gam.thr_short,Ns(2),sweep_gam.thr_long);
fprintf('=> mean-square stable when nu_bar > %.3f AND gamma_bar > %.3f.\n', crit, crit);

q={'analytic_1_minus_1overA2';'det_nu_threshold';'det_gamma_threshold'; ...
   'mc_nu_threshold_N100';'mc_nu_threshold_N300';'mc_gamma_threshold_N100';'mc_gamma_threshold_N300'};
v=[crit; nu_crit_det; gam_crit_det; sweep_nu.thr_short; sweep_nu.thr_long; ...
   sweep_gam.thr_short; sweep_gam.thr_long];
writetable(table(q,v,'VariableNames',{'quantity','value'}), ...
           fullfile(figDir,'exercise2_thresholds.csv'));

fprintf('\nExercise 2 complete. Figures and exercise2_thresholds.csv saved in: %s\n', figDir);


%% ========================= local functions =============================

function out = simulate_tcp_case(p)
% One realisation of the scalar TCP-like LQG loop (see Appendix A header).
    A=p.A; B=p.B; C=p.C; N=p.N; Qxx=p.Qxx; Quu=p.Quu;
    if isfield(p,'L') && ~isempty(p.L)
        L=p.L;
        if isfield(p,'gain_converged'); gc=p.gain_converged; else; gc=NaN; end
    else
        [L,~,gc]=tcp_control_gain(p);
    end
    nu=double(rand(1,N)<p.nu_bar); gamma=double(rand(1,N)<p.gamma_bar);
    w=sqrt(p.Q)*randn(1,N); z=sqrt(p.R)*randn(1,N);
    if isempty(p.x0); x0=sqrt(p.P0)*randn; else; x0=p.x0; end
    x=zeros(1,N+1); x(1)=x0; xhat=zeros(1,N+1); u=zeros(1,N); y=zeros(1,N);
    Ppost=zeros(1,N+1); Ppost(1)=p.P0; mu0=0;
    y(1)=gamma(1)*C*x(1)+z(1); s0=C*p.P0*C'+p.R;
    if gamma(1)==1 && s0>0
        K0=(p.P0*C')/s0; xhat(1)=mu0+K0*(y(1)-C*mu0); Ppost(1)=p.P0-K0*C*p.P0;
    else
        xhat(1)=mu0; Ppost(1)=p.P0;
    end
    u(1)=-L*xhat(1); x(2)=A*x(1)+nu(1)*B*u(1)+w(1);
    for k=1:N-1
        j=k+1; y(j)=gamma(j)*C*x(j)+z(j);
        [xhat(j),Ppost(j)]=kalman_tcp_update(p,y(j),xhat(k),u(k),nu(k),gamma(j),Ppost(k));
        u(j)=-L*xhat(j); x(j+1)=A*x(j)+nu(j)*B*u(j)+w(j);
    end
    ee=x(1:N)-xhat(1:N);
    stage=Qxx*x(1:N).^2+nu.*(Quu*u.^2);
    J_realised=Qxx*x(N+1)^2+sum(stage);
    out.x=x; out.ee=ee; out.nu=nu; out.gamma=gamma; out.gain_converged=gc;
    out.J_realised=J_realised; out.mean_cost=J_realised/N;
    out.emp_nu=mean(nu); out.emp_gamma=mean(gamma);
    out.rms_x=sqrt(mean(x.^2)); out.rms_e=sqrt(mean(ee.^2));
    out.max_absx=max(abs(x)); out.xN2=x(N+1)^2; out.P_final=Ppost(N);
end

function [xhat_upd,P_upd] = kalman_tcp_update(p,y,xhat_prev,u_prev,nu_prev,gamma_cur,P_prev)
% TCP-like Kalman update: prediction uses the known previous actuation outcome
% nu_prev; the current sensing arrival gamma_cur gates the correction.
    A=p.A; B=p.B; C=p.C; Q=p.Q; R=p.R;
    xhat_pred=A*xhat_prev+nu_prev*B*u_prev; P_pred=A*P_prev*A'+Q;
    innov=C*P_pred*C'+R;
    if gamma_cur==1 && all(innov(:)>0)
        K=(P_pred*C')/innov; xhat_upd=xhat_pred+K*(y-C*xhat_pred);
        P_upd=P_pred-K*C*P_pred; P_upd=0.5*(P_upd+P_upd.');
    else
        xhat_upd=xhat_pred; P_upd=P_pred;
    end
end

function [L,S,converged] = tcp_control_gain(p)
% Steady-state TCP-like LQG gain from the control Riccati recursion
%   L=(B'SB+Quu)^{-1}B'SA ,  S=A'SA+Qxx-nu_bar A'SB L .
    A=p.A; B=p.B; Qxx=p.Qxx; Quu=p.Quu; nu=p.nu_bar;
    S=Qxx; L=0; converged=false; Scap=1e12;
    for it=1:2000
        Ln=(B'*S*B+Quu)\(B'*S*A); Sn=A'*S*A+Qxx-nu*(A'*S*B*Ln);
        if abs(Ln-L)<=1e-12*(1+abs(Ln)); L=Ln; S=Sn; converged=true; return; end
        L=Ln; S=Sn;
        if ~isfinite(S)||S>Scap; converged=false; return; end
    end
end

function [Sss,bounded] = ctrl_steady(p)
% Steady-state control value S, iterated long enough to resolve the threshold.
% bounded is false when nu_bar is below the critical arrival probability.
    A=p.A; B=p.B; Qxx=p.Qxx; Quu=p.Quu; nu=p.nu_bar; S=Qxx; cap=1e8; bounded=true;
    for it=1:1e5
        Ln=(B'*S*B+Quu)\(B'*S*A); Sn=A'*S*A+Qxx-nu*(A'*S*B*Ln);
        if abs(Sn-S)<=1e-12*(1+abs(Sn)); Sss=Sn; bounded=true; return; end
        S=Sn;
        if ~isfinite(S)||S>cap; Sss=S; bounded=false; return; end
    end
    Sss=S; bounded=isfinite(S)&&S<=cap;
end

function tf = ctrl_bounded(p)
% True iff a stabilising steady-state controller exists (control recursion bounded).
    [~,tf]=ctrl_steady(p);
end

function [Pbar,bounded] = tcp_expected_covariance(p,maxiter,cap)
% Steady-state expected error-covariance recursion (intermittent observations):
%   Pbar <- A Pbar A' + Q - gamma_bar (A Pbar C')(C Pbar C'+R)^{-1}(C Pbar A').
% Finite fixed point iff gamma_bar > 1 - 1/A^2.
    if nargin<2||isempty(maxiter); maxiter=1e5; end
    if nargin<3||isempty(cap); cap=1e8; end
    A=p.A; C=p.C; Q=p.Q; R=p.R; g=p.gamma_bar; Pp=Q; bounded=true;
    for it=1:maxiter
        innov=C*Pp*C'+R; corr=g*(A*Pp*C')*(innov\(C*Pp*A'));
        Pn=A*Pp*A'+Q-corr;
        if abs(Pn-Pp)<=1e-12*(1+abs(Pn)); Pbar=Pn; bounded=true; return; end
        Pp=Pn;
        if ~isfinite(Pp)||Pp>cap; Pbar=Pp; bounded=false; return; end
    end
    Pbar=Pp; bounded=isfinite(Pp)&&Pp<=cap;
end

function tf = cov_bounded(p)
% True iff the expected error-covariance recursion stays bounded.
    [~,tf]=tcp_expected_covariance(p,1e5,1e8);
end

function thr = bisect_threshold(testFun, lo, hi, tol)
% Smallest value in [lo,hi] for which the monotone testFun is true (bisection).
    if ~testFun(hi); thr=NaN; return; end
    if  testFun(lo); thr=lo;  return; end
    while hi-lo>tol; mid=0.5*(lo+hi); if testFun(mid); hi=mid; else; lo=mid; end; end
    thr=hi;
end

function s = run_sweep(P, gf, which, fixedVal, Ns, masterSeed)
% One-dimensional Monte Carlo sweep of one arrival probability at two horizons,
% returning the log mean state energy and a finite-horizon growth threshold.
    M=150; logE=nan(numel(Ns),numel(gf)); grow=nan(numel(Ns),numel(gf));
    for hi=1:numel(Ns)
        p=P; p.N=Ns(hi);
        for i=1:numel(gf)
            if strcmp(which,'nu'); p.nu_bar=gf(i); p.gamma_bar=fixedVal;
            else; p.gamma_bar=gf(i); p.nu_bar=fixedVal; end
            [p.L,~,p.gain_converged]=tcp_control_gain(p);
            sumX2=zeros(1,p.N+1);
            for t=1:M
                rng(masterSeed+9000+311*hi+17*i+t);
                o=simulate_tcp_case(p); sumX2=sumX2+o.x.^2;
            end
            V=sumX2/M; logE(hi,i)=log10(max(mean(V),realmin));
            midI=round(p.N/2)+1; grow(hi,i)=V(end)/max(V(midI),realmin);
        end
    end
    s.gf=gf; s.Ns=Ns; s.logE=logE; s.grow=grow;
    s.thr_short=first_stable(gf,grow(1,:)<3); s.thr_long=first_stable(gf,grow(2,:)<3);
    s.which=which;
end

function thr = first_stable(grid, stableRow)
% Smallest grid value at which the row becomes (and stays) stable.
    thr=NaN;
    for i=1:numel(grid)
        if stableRow(i) && all(stableRow(i:end)); thr=grid(i); return; end
    end
end

function draw_heatmap(nb, gb, Z, crit, cbarLabel, useColorbar)
% Heatmap over (nu_bar, gamma_bar) with the analytic critical lines 1 - 1/A^2.
    imagesc(nb,gb,Z); set(gca,'YDir','normal'); hold on;
    xlabel('\nu_{bar} (actuation arrival)'); ylabel('\gamma_{bar} (sensing arrival)');
    if useColorbar; cb=colorbar; cb.Label.String=cbarLabel; end
    plot([crit crit],[gb(1) gb(end)],'w--','LineWidth',1.5);
    plot([nb(1) nb(end)],[crit crit],'w--','LineWidth',1.5);
    axis([nb(1) nb(end) gb(1) gb(end)]);
end

function plot_sweep(s, crit, xlab)
% log mean state energy vs arrival probability for both horizons.
    plot(s.gf,s.logE(1,:),'LineWidth',1.5); hold on; plot(s.gf,s.logE(2,:),'LineWidth',1.5);
    yl=ylim; plot([crit crit],yl,'k--','LineWidth',1.2); grid on;
    xlabel(xlab); ylabel('log_{10} mean state energy');
    legend({sprintf('N=%d',s.Ns(1)),sprintf('N=%d',s.Ns(2)),sprintf('1-1/A^2=%.3f',crit)}, ...
           'Location','northeast');
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

function C = seqmap()
% Sequential light-blue to deep-navy colour map (perceptually increasing and
% still legible when printed in grey scale). Used for the energy landscape.
    a=[0.984 0.992 1.000; 0.74 0.85 0.93; 0.33 0.58 0.79; 0.13 0.34 0.55; 0.03 0.16 0.33];
    C=interp1(linspace(0,1,size(a,1)), a, linspace(0,1,256));
end

function p = pct(x,q)
% Percentile of x at q in [0,100], computed by sorting (no toolbox needed).
    x=sort(x(:)); n=numel(x); p=x(max(1,min(n,round(q/100*n))));
end

function Zs = smooth2(Z)
% Edge-correct 3x3 moving average for display smoothing of the energy grid.
    K=ones(3); Zs=conv2(Z,K,'same')./conv2(ones(size(Z)),K,'same');
end
