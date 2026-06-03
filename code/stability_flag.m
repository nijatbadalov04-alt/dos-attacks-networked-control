function flag = stability_flag(termEnergy, rmsState)
%STABILITY_FLAG  Heuristic finite-horizon stability label from Monte Carlo
%   summaries of a case.
%
%       termEnergy : a representative (e.g. median over trials) terminal state
%                    energy x_N^2
%       rmsState   : mean RMS state over the horizon (MC mean)
%
%   Rationale. A mean-square-stable closed loop settles to an O(1-10) state
%   energy, whereas an unstable mode inflates x_N^2 by many orders of magnitude
%   over the horizon. The cutoffs below are deliberately generous and are
%   stated in the report; they should always be sanity-checked against the
%   per-case plots, because a finite horizon can only APPROXIMATE the
%   asymptotic (N -> inf) stability boundary.

    if ~isfinite(termEnergy) || termEnergy > 1e3 || rmsState > 1e2
        flag = 'unstable';
    elseif termEnergy > 25 || rmsState > 10
        flag = 'marginal';
    else
        flag = 'stable';
    end
end
