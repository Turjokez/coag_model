% run_mass_balance_biovolume_disagg_compare.m
% Compare biovolume mass-balance residuals with disaggregation OFF vs ON.

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

base.sinking_law = 'kriest_8';
base.box_depth = 2000;

% Disaggregation parameters (used when enabled)
base.c3 = 0.02;
base.r_to_rg = 1.6;
base.c4 = 1.45;

% Optional: plot rate terms by size (slower, more figures)
do_rate_terms = false;
times = [0 5 10 20 30 40 50];

% ---- output dir ----
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
outdir = fullfile(proj_root, 'output', 'figures');

% ---- run OFF ----
cfg_off = copy(base);
cfg_off.enable_disagg = false;
tag_off = sprintf('%s_disagg_off_c3_%.4g', cfg_off.sinking_law, cfg_off.c3);
[budget_off, sim_off, res_off] = run_case(cfg_off, outdir, tag_off, true, do_rate_terms, times);

% ---- run ON ----
cfg_on = copy(base);
cfg_on.enable_disagg = true;
tag_on = sprintf('%s_disagg_on_c3_%.4g', cfg_on.sinking_law, cfg_on.c3);
[budget_on, sim_on, res_on] = run_case(cfg_on, outdir, tag_on, true, do_rate_terms, times);

% ---- comparison plot ----
fig = figure('Name', 'Mass balance residual comparison');
if isfield(budget_off, 'residual_int') && ~isempty(budget_off.residual_int)
    plot(budget_off.t, budget_off.residual_int, 'LineWidth', 1.5);
    hold on;
    plot(budget_on.t, budget_on.residual_int, 'LineWidth', 1.5);
    ylabel('Biovolume (state units)');
    title(sprintf('Integrated mass-balance residual (M - M_{pred}), c3=%.4g', cfg_on.c3));
else
    plot(budget_off.t, budget_off.residual, 'LineWidth', 1.5);
    hold on;
    plot(budget_on.t, budget_on.residual, 'LineWidth', 1.5);
    ylabel('Biovolume rate (state units/day)');
    title(sprintf('Mass-balance residual (dM/dt - sum RHS), c3=%.4g', cfg_on.c3));
end
grid on;
xlabel('Time (days)');
legend('Disagg OFF', 'Disagg ON', 'Location', 'best');

if ~exist(outdir, 'dir')
    mkdir(outdir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
saveas(fig, fullfile(outdir, sprintf('mass_balance_residual_compare_c3_%.4g_%s.png', cfg_on.c3, timestamp)));

% ---- print residual summaries ----
print_summary('Disagg OFF', budget_off);
print_summary('Disagg ON', budget_on);

% ---- delta_t sensitivity (OFF only) ----
delta_ts = [1.0, 0.25];
budgets_dt = cell(size(delta_ts));
labels_dt = cell(size(delta_ts));

for i = 1:numel(delta_ts)
    dt = delta_ts(i);
    label = sprintf('Disagg OFF (dt=%.2g)', dt);

    if abs(dt - base.delta_t) < 1e-12
        budget_dt = budget_off;
    else
        cfg_dt = copy(base);
        cfg_dt.delta_t = dt;
        cfg_dt.enable_disagg = false;
        tag_dt = sprintf('%s_disagg_off_dt_%g', cfg_dt.sinking_law, dt);
        [budget_dt, ~, ~] = run_case(cfg_dt, outdir, tag_dt, false, false, times);
    end

    budgets_dt{i} = budget_dt;
    labels_dt{i} = label;
    print_summary(label, budget_dt);
end

% ---- delta_t comparison plot (OFF only) ----
use_integrated = true;
for i = 1:numel(budgets_dt)
    if ~isfield(budgets_dt{i}, 'residual_int') || isempty(budgets_dt{i}.residual_int)
        use_integrated = false;
        break;
    end
end

fig_dt = figure('Name', 'Mass balance residual vs dt (Disagg OFF)');
for i = 1:numel(budgets_dt)
    b = budgets_dt{i};
    if use_integrated
        plot(b.t, b.residual_int, 'LineWidth', 1.5);
        ylabel('Biovolume (state units)');
        title(sprintf('Mass-balance residual vs \\Delta t (Disagg OFF): integrated, c3=%.4g', cfg_on.c3));
    else
        plot(b.t, b.residual, 'LineWidth', 1.5);
        ylabel('Biovolume rate (state units/day)');
        title(sprintf('Mass-balance residual vs \\Delta t (Disagg OFF), c3=%.4g', cfg_on.c3));
    end
    hold on;
end
grid on;
xlabel('Time (days)');
legend(labels_dt, 'Location', 'best');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
saveas(fig_dt, fullfile(outdir, sprintf('mass_balance_residual_dt_sensitivity_off_c3_%.4g_%s.png', cfg_on.c3, timestamp)));

% ---- local helpers ----
function [budget, sim, res] = run_case(cfg, outdir, tag, do_plots, do_rate_terms, times)
    sim = CoagulationSimulation(cfg);
    res = sim.run();
    budget = MassBalanceBiovolume.compute(sim, res, 'state_is_biovolume', true);
    if do_plots
        MassBalanceBiovolume.plotSummary(budget, outdir, tag);
        if do_rate_terms
            MassBalanceBiovolume.plotRateTermsBySize(sim, res, times, outdir, tag);
        end
    end
end

function print_summary(label, budget)
    res_rms = sqrt(mean(budget.residual.^2));
    res_max = max(abs(budget.residual));
    fprintf('%s: RMS=%.3e, max=%.3e\n', label, res_rms, res_max);

    if isfield(budget, 'residual_int') && ~isempty(budget.residual_int)
        res_int_max = max(abs(budget.residual_int));
        denom = max(abs(budget.inventory));
        if denom > 0
            fprintf('%s: Integrated max=%.3e (rel %.3e)\n', ...
                label, res_int_max, res_int_max/denom);
        else
            fprintf('%s: Integrated max=%.3e\n', label, res_int_max);
        end
    end

    if isfield(budget, 'rate_disagg') && ~isempty(budget.rate_disagg)
        rate_disagg = budget.rate_disagg;
        disagg_rms = sqrt(mean(rate_disagg.^2));
        disagg_max = max(abs(rate_disagg));
        disagg_int_max = NaN;
        if isfield(budget, 't') && ~isempty(budget.t)
            disagg_int = cumtrapz(budget.t, rate_disagg);
            disagg_int_max = max(abs(disagg_int));
        end
        denom = max(abs(budget.inventory));
        if denom > 0 && isfinite(disagg_int_max)
            fprintf('%s: net disagg (loss term) RMS=%.3e, max=%.3e | Integrated max=%.3e (rel %.3e)\n', ...
                label, disagg_rms, disagg_max, disagg_int_max, disagg_int_max/denom);
        elseif isfinite(disagg_int_max)
            fprintf('%s: net disagg (loss term) RMS=%.3e, max=%.3e | Integrated max=%.3e\n', ...
                label, disagg_rms, disagg_max, disagg_int_max);
        else
            fprintf('%s: net disagg (loss term) RMS=%.3e, max=%.3e\n', label, disagg_rms, disagg_max);
        end
    end
end
