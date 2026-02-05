% run_disagg_c3_sweep_full_model.m
% Sweep disaggregation strength (c3) in the full model and report stability.

clear; close all; clc;
setup_paths

% ---- base configuration ----
base = SimulationConfig();
base.t_final = 50;
base.delta_t = 0.25;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;

base.enable_coag = true;
base.enable_sinking = true;
base.enable_linear = true;
base.growth = 0;
base.enable_disagg = true;
base.sinking_law = 'kriest_8';
base.box_depth = 2000;
base.c4 = 1.45;

% Sweep values (edit as needed)
c3_values = [0, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.1, 0.2];

% Optional plot
do_plot = true;

% ---- run sweep ----
n = numel(c3_values);
resid_rel = nan(1, n);
disagg_rel = nan(1, n);
resid_max = nan(1, n);
disagg_max = nan(1, n);

fprintf('Disagg c3 sweep (full model): t_final=%g days, delta_t=%g\n', ...
    base.t_final, base.delta_t);
fprintf('c3\tResid rel\tResid max\tDisagg rel\tDisagg max\n');

for i = 1:n
    cfg = copy(base);
    cfg.c3 = c3_values(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run();
    budget = MassBalanceBiovolume.compute(sim, res, 'state_is_biovolume', true);

    denom = max(abs(budget.inventory));
    if denom <= 0
        denom = 1;
    end

    % Integrated residual relative to inventory
    if isfield(budget, 'residual_int') && ~isempty(budget.residual_int)
        resid_rel(i) = max(abs(budget.residual_int)) / denom;
    else
        resid_rel(i) = max(abs(budget.residual)) / denom;
    end
    resid_max(i) = max(abs(budget.residual));

    % Net disagg (loss term) relative to inventory
    disagg_int = trapz(budget.t, budget.rate_disagg);
    disagg_rel(i) = disagg_int / denom;
    disagg_max(i) = max(abs(budget.rate_disagg));

    fprintf('%.4g\t%.3e\t%.3e\t%.3e\t%.3e\n', ...
        cfg.c3, resid_rel(i), resid_max(i), disagg_rel(i), disagg_max(i));
end

% ---- plot ----
if do_plot
    fig = figure('Name', 'Disagg c3 sweep (full model)');
    yyaxis left;
    plot(c3_values, resid_rel, 'o-', 'LineWidth', 1.5);
    ylabel('Integrated residual (rel)');
    yyaxis right;
    plot(c3_values, disagg_rel, 's--', 'LineWidth', 1.5);
    ylabel('Net disagg (rel)');
    grid on;
    xlabel('c3');
    title('Full model: residual vs net disagg across c3');
    legend('Residual rel', 'Net disagg rel', 'Location', 'best');

    script_dir = fileparts(mfilename('fullpath'));
    proj_root = fileparts(script_dir);
    outdir = fullfile(proj_root, 'output', 'figures');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    saveas(fig, fullfile(outdir, sprintf('disagg_c3_sweep_full_%s.png', timestamp)));
end
