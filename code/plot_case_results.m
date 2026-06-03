function fig = plot_case_results(out, caseInfo, outDir)
%PLOT_CASE_RESULTS  Standard four-panel diagnostic figure for one case:
%       (1) log10 running realised cost vs k
%       (2) state x_k and estimate xhat_k vs k
%       (3) estimation error x_k - xhat_k vs k, with the +/- sqrt(P_k) band
%       (4) packet-arrival sequences gamma_k (left axis) and nu_k (right axis)
%
%   caseInfo : struct with .name (file stem, e.g. 'case1') and .title (string).
%   outDir   : folder to save the PNG/FIG into; pass '' to skip saving.

    N = out.p.N;
    fig = figure('Name', caseInfo.name, 'Color', 'w', ...
                 'Position', [80 60 780 760]);

    % (1) running cost on a log10 axis (J is monotone increasing, >= 0) -------
    subplot(4,1,1);
    plot(out.k_state, log10(max(out.J_run, realmin)), 'LineWidth', 1.3);
    grid on; xlim([0 N]);
    xlabel('k'); ylabel('log_{10} J(k)');
    title('Running realised cost');

    % (2) true state and TCP-like estimate -----------------------------------
    subplot(4,1,2);
    stairs(out.k_state, out.x, 'LineWidth', 1.1); hold on;
    stairs(out.k_io,   out.xhat, '--', 'LineWidth', 1.1);
    grid on; xlim([0 N]);
    xlabel('k'); ylabel('state');
    legend({'x_k','xhat_k'}, 'Location', 'best');
    title('State and TCP-like estimate');

    % (3) estimation error with covariance band ------------------------------
    subplot(4,1,3);
    band = sqrt(max(out.Ppost, 0));
    plot(out.k_io, out.ee, 'LineWidth', 1.1); hold on;
    plot(out.k_io,  band, ':', 'LineWidth', 1.0);
    plot(out.k_io, -band, ':', 'LineWidth', 1.0);
    grid on; xlim([0 N]);
    xlabel('k'); ylabel('x_k - xhat_k');
    title('Estimation error (dotted: \pm one std. dev. \surd P_k)', 'Interpreter','tex');

    % (4) packet arrivals ----------------------------------------------------
    subplot(4,1,4);
    yyaxis left;  stairs(out.k_io, out.gamma, 'LineWidth', 1.0);
    ylim([-0.2 1.2]); ylabel('\gamma_k');
    yyaxis right; stairs(out.k_io, out.nu, 'LineWidth', 1.0);
    ylim([-0.2 1.2]); ylabel('\nu_k');
    xlim([0 N]); xlabel('k');
    title('Packet arrivals (\gamma_k sensing, \nu_k actuation)');

    % overall title (sgtitle is R2018b+; fall back to an annotation) ---------
    if isfield(caseInfo, 'title')
        try
            sgtitle(caseInfo.title, 'FontWeight', 'bold');
        catch
            annotation(fig, 'textbox', [0 0.965 1 0.03], 'String', caseInfo.title, ...
                'HorizontalAlignment','center', 'EdgeColor','none', 'FontWeight','bold');
        end
    end

    if ~isempty(outDir)
        save_figure(fig, outDir, caseInfo.name);
    end
end
