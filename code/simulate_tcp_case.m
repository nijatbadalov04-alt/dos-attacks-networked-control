function out = simulate_tcp_case(p)
%SIMULATE_TCP_CASE  Simulate ONE realisation of the scalar TCP-like LQG system
%   with Bernoulli actuation (nu) and sensing (gamma) packet drops.
%
%       x_{k+1} = A x_k + nu_k B u_k + w_k ,   w_k ~ N(0, Q)
%       y_k     = gamma_k C x_k + z_k       ,   z_k ~ N(0, R)
%       u_k     = -L xhat_k                 ,   xhat from the TCP-like filter
%
%   INDEXING CONVENTION  (math time k = 0..N  <->  MATLAB index k+1)
%       x(k+1)    = x_k       for k = 0..N      (N+1 states)
%       xhat(k+1) = xhat_k    for k = 0..N-1    (posterior estimates)
%       u(k+1)    = u_k       for k = 0..N-1
%       y(k+1)    = y_k       for k = 0..N-1
%       nu(k+1) / gamma(k+1) / w(k+1) / z(k+1)  hold the k-th realisation,
%                                               k = 0..N-1
%
%   All randomness uses base MATLAB only (randn / rand), so NO toolboxes are
%   required (the shared template used ss, mvnrnd and binornd). The CALLER is
%   responsible for seeding (rng) before each call, so Monte Carlo trials are
%   reproducible; this function does not touch the global rng state.
%
%   INPUT   p : struct from base_params() with A,N,Q,R,nu_bar,gamma_bar set.
%               Optional: p.L (precomputed steady-state gain) and
%               p.gain_converged, to avoid recomputing the gain in MC loops.
%   OUTPUT  out : struct with the full trajectories, packet sequences, running
%               cost and a set of scalar performance metrics (listed at end).

    % ---------- light input validation (cheap; toggled by p.validate) -------
    if ~isfield(p,'validate') || p.validate
        assert(p.nu_bar    >= 0 && p.nu_bar    <= 1, 'nu_bar must be in [0,1]');
        assert(p.gamma_bar >= 0 && p.gamma_bar <= 1, 'gamma_bar must be in [0,1]');
        assert(p.N >= 2 && p.N == round(p.N),        'N must be an integer >= 2');
        assert(p.Q >= 0 && p.R >= 0,                 'noise variances must be >= 0');
    end

    A = p.A; B = p.B; C = p.C; N = p.N; Qxx = p.Qxx; Quu = p.Quu;

    % ---------- steady-state control gain (reuse if provided) ---------------
    if isfield(p,'L') && ~isempty(p.L)
        L = p.L;
        if isfield(p,'gain_converged'); gain_converged = p.gain_converged;
        else;                           gain_converged = NaN; end
    else
        [L, ~, gain_converged] = tcp_control_gain(p);
    end

    % ---------- draw all randomness up front (k = 0..N-1 -> idx 1..N) -------
    nu    = double(rand(1,N) < p.nu_bar);     % actuation arrivals  nu_k
    gamma = double(rand(1,N) < p.gamma_bar);  % sensing  arrivals   gamma_k
    w     = sqrt(p.Q)*randn(1,N);             % process noise       w_k
    z     = sqrt(p.R)*randn(1,N);             % measurement noise   z_k
    if isempty(p.x0)
        x0 = sqrt(p.P0)*randn;                % x_0 ~ N(0,P0)
    else
        x0 = p.x0;
    end

    % ---------- allocate ----------------------------------------------------
    x     = zeros(1,N+1);  x(1)     = x0;     % x_0..x_N
    xhat  = zeros(1,N+1);                     % xhat_0..xhat_{N-1} used
    u     = zeros(1,N);                       % u_0..u_{N-1}
    y     = zeros(1,N);                       % y_0..y_{N-1}
    Ppost = zeros(1,N+1);  Ppost(1) = p.P0;   % posterior covariance per k
    Ppred = zeros(1,N+1);  Ppred(1) = p.P0;
    Kk    = zeros(1,N);

    mu0 = 0;                                  % prior mean of x_0

    % ---------- k = 0: prior, gated correction with y_0, then u_0 -----------
    y(1)   = gamma(1)*C*x(1) + z(1);
    innov0 = C*p.P0*C' + p.R;
    if gamma(1) == 1 && innov0 > 0
        K0       = (p.P0*C')/innov0;
        xhat(1)  = mu0 + K0*(y(1) - C*mu0);
        Ppost(1) = p.P0 - K0*C*p.P0;
        Kk(1)    = K0;
    else
        xhat(1)  = mu0;
        Ppost(1) = p.P0;
    end
    u(1) = -L*xhat(1);                         % u_0
    x(2) = A*x(1) + nu(1)*B*u(1) + w(1);       % x_1

    % ---------- k = 1..N-1: TCP predict/correct, control, propagate ---------
    for k = 1:N-1
        j    = k + 1;                          % array index for time k
        y(j) = gamma(j)*C*x(j) + z(j);         % measurement y_k
        [xhat(j), Ppost(j), ~, Ppred(j), Kk(j)] = kalman_tcp_update( ...
            p, y(j), xhat(k), u(k), nu(k), gamma(j), Ppost(k));
        u(j)   = -L*xhat(j);                   % control u_k
        x(j+1) = A*x(j) + nu(j)*B*u(j) + w(j); % state x_{k+1}
    end

    % ---------- estimation error e_k = x_k - xhat_k, k = 0..N-1 -------------
    ee = x(1:N) - xhat(1:N);

    % ---------- realised cost (eq. (3) without the expectation) ------------
    %   J = Qxx x_N^2 + sum_{k=0}^{N-1} ( Qxx x_k^2 + nu_k Quu u_k^2 )
    stage      = Qxx*x(1:N).^2 + nu.*(Quu*u.^2);     % per-step stage cost
    J_run      = cumsum(stage) + Qxx*x(2:N+1).^2;    % running cost J(m), m=1..N
    J_run      = [Qxx*x(1)^2, J_run];                % prepend m = 0 term
    J_realised = Qxx*x(N+1)^2 + sum(stage);

    % ---------- pack output -------------------------------------------------
    out                = struct();
    out.p              = p;
    out.k_state        = 0:N;            % time axis for states (length N+1)
    out.k_io           = 0:N-1;          % time axis for inputs/meas/errors
    out.x              = x;
    out.xhat           = xhat(1:N);      % defined for k = 0..N-1
    out.ee             = ee;
    out.u              = u;
    out.y              = y;
    out.nu             = nu;
    out.gamma          = gamma;
    out.Ppost          = Ppost(1:N);
    out.Ppred          = Ppred(1:N);
    out.K              = Kk;
    out.J_run          = J_run;          % J_run(m+1) = cost accumulated to time m
    out.L              = L;
    out.gain_converged = gain_converged;

    % scalar performance metrics
    out.J_realised  = J_realised;
    out.mean_cost   = J_realised / N;        % per-step average cost
    out.emp_nu      = mean(nu);              % realised mean actuation arrival
    out.emp_gamma   = mean(gamma);           % realised mean sensing arrival
    out.rms_x       = sqrt(mean(x.^2));      % RMS state over k = 0..N
    out.rms_e       = sqrt(mean(ee.^2));     % RMS estimation error over k=0..N-1
    out.max_absx    = max(abs(x));           % peak |state|
    out.xN2         = x(N+1)^2;              % terminal state energy (divergence probe)
    out.P_final     = Ppost(N);             % last posterior covariance computed
end
