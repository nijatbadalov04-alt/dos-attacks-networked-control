function [Pbar, bounded] = tcp_expected_covariance(p, maxiter, cap)
%TCP_EXPECTED_COVARIANCE  Steady-state EXPECTED (mean-field) prediction error
%   covariance for the TCP-like filter with i.i.d. sensing dropouts gamma_bar.
%
%   This iterates the modified algebraic Riccati equation for estimation with
%   intermittent observations (Sinopoli et al., 2004) - the same recursion as
%   the P_hat upper bound in the shared template DoS_control_TCP_26.m:
%
%       Pbar <- A Pbar A' + Q
%               - gamma_bar * (A Pbar C') (C Pbar C' + R)^{-1} (C Pbar A')
%
%   For the scalar unstable plant this recursion has a finite fixed point iff
%       gamma_bar > 1 - 1/A^2 ,
%   and Pbar diverges below that critical sensing arrival probability. The
%   recursion is DETERMINISTIC, so it locates the estimator stability boundary
%   without Monte Carlo noise - the right tool for comparing with theory.
%
%   Inputs : p with A,C,Q,R,gamma_bar; optional maxiter, cap.
%   Outputs: Pbar    steady-state expected covariance (cap-limited if diverging)
%            bounded true iff a finite fixed point was reached / stayed bounded.

    if nargin < 2 || isempty(maxiter); maxiter = 1e5; end
    if nargin < 3 || isempty(cap);     cap     = 1e8; end

    A = p.A; C = p.C; Q = p.Q; R = p.R; g = p.gamma_bar;

    P = Q;                          % start from the process-noise level
    bounded = true;
    for it = 1:maxiter
        innov = C*P*C' + R;
        corr  = g*(A*P*C') * (innov \ (C*P*A'));   % no inv(); right/left division
        Pn    = A*P*A' + Q - corr;
        if abs(Pn - P) <= 1e-12*(1 + abs(Pn))
            Pbar = Pn; bounded = true; return;     % converged fixed point
        end
        P = Pn;
        if ~isfinite(P) || P > cap
            Pbar = P; bounded = false; return;     % diverging (gamma below crit)
        end
    end
    Pbar = P; bounded = isfinite(P) && P <= cap;    % bounded after maxiter
end
