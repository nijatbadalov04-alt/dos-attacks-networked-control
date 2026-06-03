function p = base_params()
%BASE_PARAMS  Common configuration for the scalar TCP-like DoS lab system.
%
%   p = BASE_PARAMS() returns a struct with the fields shared by both
%   exercises of the ELE419 DoS laboratory. Each exercise overrides the few
%   fields that differ (A, N, the process/measurement noise variances, and the
%   packet-arrival probabilities nu_bar / gamma_bar).
%
%   Field meanings (notation follows the assignment brief, eqs. (1)-(3)):
%       A,B,C      scalar state-space coefficients in
%                      x_{k+1} = A x_k + nu_k B u_k + w_k
%                      y_k     = gamma_k C x_k + z_k
%       N          horizon: inputs/measurements for k = 0..N-1, states x_0..x_N
%       Qxx,Quu    LQ COST penalties on state and input (identity = 1 per brief)
%       Q          process-noise VARIANCE  Var(w_k)   (kept named Q as template)
%       R          measurement-noise VARIANCE Var(z_k)
%       nu_bar     P[nu_k = 1]    actuation packet-ARRIVAL probability
%       gamma_bar  P[gamma_k = 1] sensing  packet-ARRIVAL probability
%       P0         initial estimation-error covariance
%       x0         [] -> drawn ~ N(0,P0) each trial; or a fixed scalar
%       validate   true -> run cheap input checks inside simulate_tcp_case
%
%   NAMING NOTE. The brief overloads "Q": Qxx/Quu are COST matrices while the
%   noise is N(0, var). We therefore keep Q and R for the NOISE variances (as
%   in the shared template DoS_control_TCP_26.m) and use Qxx,Quu for the COST
%   penalties. This avoids any ambiguity in the code.

    p           = struct();
    p.B         = 1;       % actuator coefficient (scalar)
    p.C         = 1;       % sensor coefficient (scalar)
    p.Qxx       = 1;       % state penalty in the cost (identity per brief)
    p.Quu       = 1;       % input penalty in the cost (identity per brief)
    p.P0        = 0.01;    % small initial error covariance (matches template)
    p.x0        = [];      % [] => random initial state x_0 ~ N(0,P0)
    p.validate  = true;    % cheap input validation in the simulator
end
