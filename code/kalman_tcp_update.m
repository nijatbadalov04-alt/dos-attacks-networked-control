function [xhat_upd, P_upd, xhat_pred, P_pred, K] = ...
        kalman_tcp_update(p, y, xhat_prev, u_prev, nu_prev, gamma_cur, P_prev)
%KALMAN_TCP_UPDATE  One step of the TCP-like Kalman filter.
%
%   This is a clean, toolbox-free rewrite of the shared KalmanTCP.m. It keeps
%   exactly the same mathematics but uses matrix right-division instead of
%   inv(), enforces covariance symmetry, and handles the noiseless-measurement
%   case (R = 0) safely.
%
%   TCP-like information set  F_k = { y^k, gamma^k, nu^{k-1} }: the controller
%   receives an acknowledgement of the PREVIOUS actuation outcome nu_{k-1}, so
%   it enters the prediction deterministically; the CURRENT sensing arrival
%   gamma_k gates the measurement correction.
%
%       Prediction : xhat^-_k = A*xhat_{k-1} + nu_{k-1}*B*u_{k-1}
%                    P^-_k     = A*P_{k-1}*A' + Q
%       Correction : if gamma_k = 1   (measurement packet received)
%                       K_k    = P^-_k C' (C P^-_k C' + R)^{-1}
%                       xhat_k = xhat^-_k + K_k (y_k - C xhat^-_k)
%                       P_k    = P^-_k - K_k C P^-_k
%                    else            (no measurement -> propagate prediction)
%                       xhat_k = xhat^-_k ,  P_k = P^-_k
%
%   Inputs
%       p          struct with A,B,C,Q (process-noise var), R (meas-noise var)
%       y          current measurement y_k
%       xhat_prev  previous posterior estimate xhat_{k-1}
%       u_prev     previous input u_{k-1}
%       nu_prev    previous actuation outcome nu_{k-1}  (KNOWN via TCP ack)
%       gamma_cur  current sensing arrival gamma_k (0/1)
%       P_prev     previous posterior covariance P_{k-1}
%
%   WHY R = 0 IS SAFE. The innovation covariance is S = C P^-_k C' + R. Since
%   P^-_k = A P_{k-1} A' + Q >= Q > 0 throughout this lab, S stays strictly
%   positive even when R = 0, so no singular inverse occurs and the noiseless
%   case simply yields K = 1 (scalar), i.e. xhat_k = y_k = x_k. A defensive
%   guard still skips the correction if S were ever non-positive.

    A = p.A; B = p.B; C = p.C; Q = p.Q; R = p.R;

    % --- Prediction (uses the known previous actuation outcome nu_prev) ---
    xhat_pred = A*xhat_prev + nu_prev*B*u_prev;
    P_pred    = A*P_prev*A' + Q;

    % --- Correction (only when the measurement packet arrives) ---
    innov_cov = C*P_pred*C' + R;                 % innovation covariance S
    if gamma_cur == 1 && all(innov_cov(:) > 0)
        K        = (P_pred*C') / innov_cov;      % Kalman gain (no inv())
        xhat_upd = xhat_pred + K*(y - C*xhat_pred);
        P_upd    = P_pred - K*C*P_pred;
        P_upd    = 0.5*(P_upd + P_upd.');        % symmetrise (matrix-safe)
    else
        K        = zeros(size(P_pred*C'));
        xhat_upd = xhat_pred;                    % predict-only
        P_upd    = P_pred;
    end
end
