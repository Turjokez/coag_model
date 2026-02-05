% run_disagg_c3_sweep_conservation.m
% Sweep disaggregation strength (c3) and report inventory drift.

clear; close all; clc;
setup_paths

% ---- base configuration ----
base = SimulationConfig();
base.t_final = 50;
base.delta_t = 0.25;
base.n_sections = 35;

% Disaggregation-only setup
base.enable_pp = false;
base.pp_source = 0;
base.enable_coag = false;
base.enable_sinking = false;
base.enable_linear = false;
base.growth = 0;
base.enable_disagg = true;
base.c4 = 1.45;

% Sweep values (edit as needed)
c3_values = [0, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05];
target_drift_rel = -0.05; % target relative inventory drift (e.g., -5%)

% Optional plot
do_plot = true;

% ---- run sweep ----
drift_rel = zeros(size(c3_values));
disagg_int_rel = zeros(size(c3_values));

fprintf('Disagg c3 sweep (disagg-only): t_final=%g days, delta_t=%g\n', ...
    base.t_final, base.delta_t);
fprintf('c3\tDrift rel\tNet disagg rel\n');

for i = 1:numel(c3_values)
    cfg = copy(base);
    cfg.c3 = c3_values(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run();
    budget = MassBalanceBiovolume.compute(sim, res, 'state_is_biovolume', true);

    M = budget.inventory;
    denom = max(abs(M));
    if denom <= 0
        denom = 1;
    end

    drift = M(end) - M(1);
    drift_rel(i) = drift / denom;

    disagg_int = trapz(budget.t, budget.rate_disagg);
    disagg_int_rel(i) = disagg_int / denom;

    fprintf('%.3g\t%.3e\t%.3e\n', cfg.c3, drift_rel(i), disagg_int_rel(i));
end

% ---- report closest to target ----
[~, idx] = min(abs(drift_rel - target_drift_rel));
fprintf('Target drift %.3f: closest c3=%.3g (drift rel %.3e)\n', ...
    target_drift_rel, c3_values(idx), drift_rel(idx));

% ---- plot ----
if do_plot
    fig = figure('Name', 'Disagg c3 sweep: inventory drift');
    plot(c3_values, drift_rel, 'o-', 'LineWidth', 1.5);
    grid on;
    xlabel('c3');
    ylabel('Inventory drift (rel)');
    title('Disaggregation-only drift vs c3');

    script_dir = fileparts(mfilename('fullpath'));
    proj_root = fileparts(script_dir);
    outdir = fullfile(proj_root, 'output', 'figures');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    saveas(fig, fullfile(outdir, sprintf('disagg_c3_sweep_drift_%s.png', timestamp)));
end
