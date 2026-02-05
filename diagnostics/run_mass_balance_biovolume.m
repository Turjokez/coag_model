% run_mass_balance_biovolume.m
% Biovolume-consistent mass-balance diagnostics (Adrian-style)

clear; close all; clc;
setup_paths

% ---- configuration (edit if needed) ----
cfg = SimulationConfig();
cfg.t_final = 50;
cfg.delta_t = 0.25;
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = false;
cfg.enable_linear = true;
cfg.growth = 0;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

% ---- run ----
sim = CoagulationSimulation(cfg);
res = sim.run();

% ---- compute budget ----
% State vector in this model is already biovolume/volume units
budget = MassBalanceBiovolume.compute(sim, res, 'state_is_biovolume', true);

% ---- plots ----
% Resolve project root (one level above diagnostics/)
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
outdir = fullfile(proj_root, 'output', 'figures');
MassBalanceBiovolume.plotSummary(budget, outdir, 'kriest8');

% Rate terms by size at selected times (edit as needed)
times = [0 5 10 20 30 40 50];
MassBalanceBiovolume.plotRateTermsBySize(sim, res, times, outdir, 'kriest8');

% ---- print residual summary ----
res_rms = sqrt(mean(budget.residual.^2));
res_max = max(abs(budget.residual));
fprintf('Biovolume budget residual: RMS=%.3e, max=%.3e\n', res_rms, res_max);

% Integrated residual summary (more robust than finite-diff)
if isfield(budget, 'residual_int')
    res_int_max = max(abs(budget.residual_int));
    denom = max(abs(budget.inventory));
    if denom > 0
        fprintf('Integrated residual: max=%.3e (rel %.3e)\n', res_int_max, res_int_max/denom);
    else
        fprintf('Integrated residual: max=%.3e\n', res_int_max);
    end
end

% ---- auto-save all open figures (PNG) ----
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
tag = 'kriest8';
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
figs = findall(0, 'Type', 'figure');
for k = 1:numel(figs)
    fig = figs(k);
    fname = get(fig, 'Name');
    if isempty(fname)
        fname = sprintf('figure_%d', fig.Number);
    end
    fname = regexprep(fname, '[^A-Za-z0-9_\\-]', '_');

    % Hide toolbar/menubar to avoid export toolbar in images
    old_tb = get(fig, 'ToolBar');
    old_mb = get(fig, 'MenuBar');
    set(fig, 'ToolBar', 'none');
    set(fig, 'MenuBar', 'none');

    outpath = fullfile(outdir, sprintf('%s_%s_%s.png', fname, tag, timestamp));
    saveas(fig, outpath);

    % Restore UI
    set(fig, 'ToolBar', old_tb);
    set(fig, 'MenuBar', old_mb);
end
fprintf('Saved %d figure(s) to %s\n', numel(figs), outdir);
