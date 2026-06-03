function save_figure(figHandle, outDir, baseName)
%SAVE_FIGURE  Save a figure as PNG (always) and FIG (best effort), with
%   fallbacks so the project works across MATLAB versions and even without the
%   newer graphics functions.
%
%   save_figure(figHandle, outDir, baseName) writes
%       <outDir>/<baseName>.png   and   <outDir>/<baseName>.fig

    if ~exist(outDir, 'dir'); mkdir(outDir); end
    pngPath = fullfile(outDir, [baseName '.png']);
    figPath = fullfile(outDir, [baseName '.fig']);

    % --- print-friendly appearance ---
    %   Force a LIGHT (white) theme so figures print well, regardless of the
    %   MATLAB session's default light/dark theme (R2025a+ uses themes), and
    %   hide the per-axes toolbar so exportgraphics does not warn about it.
    try; theme(figHandle, 'light'); catch; end       %#ok<*CTCH>  (older releases)
    set(figHandle, 'Color', 'w');
    axList = findall(figHandle, 'Type', 'axes');
    for a = reshape(axList, 1, [])
        try; a.Toolbar.Visible = 'off'; catch; end
    end

    % --- PNG: prefer exportgraphics (R2020a+), fall back to print/saveas ---
    try
        exportgraphics(figHandle, pngPath, 'Resolution', 150);
    catch
        try
            print(figHandle, pngPath, '-dpng', '-r150');
        catch
            saveas(figHandle, pngPath);
        end
    end

    % --- FIG: best effort (skip silently if unsupported) ---
    try
        savefig(figHandle, figPath);
    catch
        % no .fig produced; PNG is still written
    end
end
