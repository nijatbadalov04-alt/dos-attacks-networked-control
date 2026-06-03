function main_ELE419_DoS_lab()
%MAIN_ELE419_DoS_LAB  One-command driver for the ELE419 DoS laboratory.
%
%   Running this file from a clean MATLAB session reproduces EVERYTHING:
%       * creates the outputs/ folder tree,
%       * runs Exercise 1 (six DoS cases, A = 1.1) with Monte Carlo support,
%       * runs Exercise 2 (stability search over nu_bar, gamma_bar, A = 1.255),
%       * saves all figures (PNG + FIG) and tables (CSV + MAT + LaTeX),
%       * writes a timestamped run log and prints concise conclusions.
%
%   No toolboxes are required: the simulation core uses only base MATLAB
%   (randn/rand). Tested on MATLAB R2025b.
%
%   Reproducibility: a single master seed is set here; every case and every
%   Monte Carlo trial derives a deterministic sub-seed from it, so re-running
%   gives identical numbers and figures.

    close all;
    thisDir = fileparts(mfilename('fullpath'));
    cd(thisDir);
    addpath(thisDir);

    % ---------------- output folder tree -----------------------------------
    outRoot         = fullfile(thisDir, 'outputs');
    outDirs.root    = outRoot;
    outDirs.ex1     = fullfile(outRoot, 'exercise1');
    outDirs.ex2     = fullfile(outRoot, 'exercise2');
    outDirs.tables  = fullfile(outRoot, 'tables');
    outDirs.figures = fullfile(outRoot, 'figures');
    fn = fieldnames(outDirs);
    for i = 1:numel(fn)
        if ~exist(outDirs.(fn{i}), 'dir'); mkdir(outDirs.(fn{i})); end
    end

    % ---------------- run log (diary) --------------------------------------
    logFile = fullfile(outRoot, sprintf('run_log_%s.txt', ...
                       datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
    diary(logFile); diary on;
    cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>
    fprintf('ELE419 DoS lab - run started %s\n', datestr(now)); %#ok<TNOW1,DATST>
    fprintf('Output directory: %s\n', outRoot);

    % ---------------- global configuration / reproducibility ---------------
    cfg.masterSeed = 2026;          % master seed (shared template used rng(2026))
    rng(cfg.masterSeed);

    cfg.mcTrials   = 300;           % Exercise 1: Monte Carlo trials per case

    cfg.ex2grid    = 0:0.05:1;      % Exercise 2: (nu_bar, gamma_bar) grid (21x21)
    cfg.ex2trials  = 80;            % Exercise 2: MC trials per grid point
    cfg.ex2horizon = 100;           % Exercise 2: validation horizon (brief N=100)

    fprintf(['Config: masterSeed=%d, mcTrials=%d, ex2grid=%g:%g:%g, ' ...
             'ex2trials=%d, ex2horizon=%d\n'], cfg.masterSeed, cfg.mcTrials, ...
             cfg.ex2grid(1), cfg.ex2grid(2)-cfg.ex2grid(1), cfg.ex2grid(end), ...
             cfg.ex2trials, cfg.ex2horizon);

    % ---------------- run the two exercises --------------------------------
    tAll     = tic;
    results1 = run_exercise1(cfg, outDirs);
    results2 = run_exercise2(cfg, outDirs);

    save(fullfile(outDirs.tables, 'all_results.mat'), 'results1', 'results2', 'cfg');
    fprintf('\nAll done in %.1f s. Figures and tables are under:\n  %s\n', ...
            toc(tAll), outRoot);
    diary off;
end
