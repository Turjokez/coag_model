% run_disagg_only_conservation.m
% Quick check: disaggregation-only net contribution to total biovolume.

clear; close all; clc;
setup_paths

% ---- configuration ----
cfg = SimulationConfig();
cfg.t_final = 50;
cfg.delta_t = 0.25;
cfg.n_sections = 35;

% Turn everything OFF except disaggregation
cfg.enable_pp = false;
cfg.pp_source = 0;
cfg.enable_coag = false;
cfg.enable_sinking = false;
cfg.enable_linear = false;
cfg.growth = 0;
cfg.enable_disagg = true;

% Disaggregation parameters (if you want to sweep, edit here)
cfg.c3 = 0.02;
cfg.c4 = 1.45;

% Optional plots
do_plots = false;

% ---- run ----
sim = CoagulationSimulation(cfg);
res = sim.run();
budget = MassBalanceBiovolume.compute(sim, res, 'state_is_biovolume', true);

% ---- summary ----
fprintf('Disagg-only conservation check (coag/linear/sinking/PP OFF)\n');

% Residuals (should reflect disagg-only behavior)
res_rms = sqrt(mean(budget.residual.^2));
res_max = max(abs(budget.residual));
fprintf('Residual: RMS=%.3e, max=%.3e\n', res_rms, res_max);

if isfield(budget, 'residual_int') && ~isempty(budget.residual_int)
    res_int_max = max(abs(budget.residual_int));
    denom = max(abs(budget.inventory));
    if denom > 0
        fprintf('Integrated residual: max=%.3e (rel %.3e)\n', res_int_max, res_int_max/denom);
    else
        fprintf('Integrated residual: max=%.3e\n', res_int_max);
    end
end

% Net disagg contribution (instantaneous + integrated)
rate_disagg = budget.rate_disagg;
disagg_rms = sqrt(mean(rate_disagg.^2));
disagg_max = max(abs(rate_disagg));
disagg_int = trapz(budget.t, rate_disagg);
denom = max(abs(budget.inventory));
if denom > 0
    fprintf('Net disagg (loss term): RMS=%.3e, max=%.3e | Integrated=%.3e (rel %.3e)\n', ...
        disagg_rms, disagg_max, disagg_int, disagg_int/denom);
else
    fprintf('Net disagg (loss term): RMS=%.3e, max=%.3e | Integrated=%.3e\n', ...
        disagg_rms, disagg_max, disagg_int);
end

% Inventory drift check
M = budget.inventory;
delta_M = M(end) - M(1);
if denom > 0
    fprintf('Inventory drift: dM=%.3e (rel %.3e)\n', delta_M, delta_M/denom);
else
    fprintf('Inventory drift: dM=%.3e\n', delta_M);
end

% ---- optional plots ----
if do_plots
    script_dir = fileparts(mfilename('fullpath'));
    proj_root = fileparts(script_dir);
    outdir = fullfile(proj_root, 'output', 'figures');
    MassBalanceBiovolume.plotSummary(budget, outdir, 'disagg_only');
end
