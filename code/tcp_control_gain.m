function [L, S, converged] = tcp_control_gain(p)
%TCP_CONTROL_GAIN  Steady-state TCP-like LQ control gain via the modified
%   algebraic Riccati recursion (MARE) used in the shared template:
%
%       L     = (B' S B + Quu)^{-1} B' S A
%       S_new = A' S A + Qxx - nu_bar * A' S B L
%
%   The shared file DoS_control_TCP_26.m runs this SAME recursion forward in
%   time, alongside the simulation. Here we iterate it to its fixed point and
%   return the gain it settles to after the initial transient. That fixed
%   point is the standard infinite-horizon / steady-state controller, and over
%   the lab horizon (N = 100) the time-varying finite-horizon gain converges to
%   it within a few steps, so applying the constant gain u_k = -L*xhat_k is
%   both cleaner and numerically equivalent to the template after the transient.
%
%   THRESHOLD BEHAVIOUR. For the scalar unstable plant, the MARE has a finite
%   fixed point only when the actuation arrival probability exceeds the
%   critical value 1 - 1/A^2. Below it, S grows without bound; the GAIN,
%   however, still converges (L -> A/B as S -> inf). We therefore cap S, return
%   the limiting gain, and set converged=false so the caller can report that no
%   stabilising steady-state controller exists for that nu_bar.
%
%   Inputs : p with fields A,B,Qxx,Quu,nu_bar.
%   Outputs: L         steady-state (or limiting) control gain
%            S         steady-state value matrix (Inf-capped if diverging)
%            converged true if the MARE reached a finite fixed point

    A = p.A; B = p.B; Qxx = p.Qxx; Quu = p.Quu; nu = p.nu_bar;

    S = Qxx;            % initialise at the terminal cost weight
    L = 0;
    converged = false;
    Scap = 1e12;        % guard against MARE divergence (nu below critical)

    for it = 1:2000
        L_new = (B'*S*B + Quu) \ (B'*S*A);
        S_new = A'*S*A + Qxx - nu*(A'*S*B*L_new);

        if abs(L_new - L) <= 1e-12*(1 + abs(L_new))
            L = L_new; S = S_new; converged = true; return;     % fixed point
        end
        L = L_new; S = S_new;

        if ~isfinite(S) || S > Scap
            % MARE diverging: the gain has effectively reached its limit (A/B).
            converged = false; return;
        end
    end
end
