% run_disagg_operator_split_compare.m
% Compare legacy vs operator-split disaggregation in 0-D.

clear; close all; clc;
setup_paths

% --------- base configuration ---------
cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 2;      % output step (days)
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_linear = true;
cfg.growth = 0;
cfg.enable_disagg = true;

cfg.c3 = 0.02;
cfg.r_to_rg = 1.6;
cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

% Operator-split parameters
use_dmax_override = true;
dmax_override_cm = 50;  % set to [] to use epsilon instead
eps_val = 1e-6;         % turbulence dissipation rate (example)
outer_dt = 1/24;        % days (1 hour)

% --------- legacy run ---------
cfg_legacy = cfg.copy();
cfg_legacy.disagg_mode = 'legacy';

sim_legacy = CoagulationSimulation(cfg_legacy);
res_legacy = sim_legacy.run();

% --------- operator-split run ---------
cfg_split = cfg.copy();
cfg_split.disagg_mode = 'operator_split';
if use_dmax_override && ~isempty(dmax_override_cm)
    cfg_split.disagg_dmax_cm = dmax_override_cm;
    cfg_split.disagg_epsilon = [];
else
    cfg_split.disagg_epsilon = eps_val;
    cfg_split.disagg_dmax_cm = [];
end
cfg_split.disagg_outer_dt = outer_dt;
cfg_split.disagg_frac_next = 2/3;
cfg_split.disagg_C = 3.0;
cfg_split.disagg_gamma = 0.15;

sim_split = CoagulationSimulation(cfg_split);
res_split = sim_split.run();

% --------- export fractions over time ---------
t = res_legacy.time(:);
Fv_L = res_legacy.output_data.fluxspec;
Fi_L = res_legacy.output_data.fluxspec_i;
Fv_S = res_split.output_data.fluxspec;
Fi_S = res_split.output_data.fluxspec_i;

d_cm_v = res_legacy.output_data.diam_v(:);
d_cm_i = res_legacy.output_data.diam_i(:);

fS_v_L = zeros(numel(t),1); fM_v_L = fS_v_L; fL_v_L = fS_v_L;
fS_i_L = fS_v_L; fM_i_L = fS_v_L; fL_i_L = fS_v_L;
fS_v_S = fS_v_L; fM_v_S = fS_v_L; fL_v_S = fS_v_L;
fS_i_S = fS_v_L; fM_i_S = fS_v_L; fL_i_S = fS_v_L;

for k = 1:numel(t)
    [fS_v_L(k), fM_v_L(k), fL_v_L(k)] = export_size_classes(d_cm_v, Fv_L(k,:).', 'volume');
    [fS_i_L(k), fM_i_L(k), fL_i_L(k)] = export_size_classes(d_cm_i, Fi_L(k,:).', 'image');
    [fS_v_S(k), fM_v_S(k), fL_v_S(k)] = export_size_classes(d_cm_v, Fv_S(k,:).', 'volume');
    [fS_i_S(k), fM_i_S(k), fL_i_S(k)] = export_size_classes(d_cm_i, Fi_S(k,:).', 'image');
end

% --------- print comparison summary ---------
targets = [50, 100, 200];
fprintf('\n=== Export fractions (legacy vs split) ===\n');
for tt = targets
    [~,idx] = min(abs(t - tt));
    fprintf('t=%.0f d (vol)  legacy: S=%.3f M=%.3f L=%.3f | split: S=%.3f M=%.3f L=%.3f\n', ...
        t(idx), fS_v_L(idx), fM_v_L(idx), fL_v_L(idx), fS_v_S(idx), fM_v_S(idx), fL_v_S(idx));
    fprintf('t=%.0f d (img)  legacy: S=%.3f M=%.3f L=%.3f | split: S=%.3f M=%.3f L=%.3f\n', ...
        t(idx), fS_i_L(idx), fM_i_L(idx), fL_i_L(idx), fS_i_S(idx), fM_i_S(idx), fL_i_S(idx));
end

% --------- plots ---------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

if use_dmax_override && ~isempty(dmax_override_cm)
    tag = sprintf('dmax_%g', dmax_override_cm);
else
    tag = sprintf('eps_%g', eps_val);
end

fig1 = figure; hold on;
plot(t, fS_i_L, 'LineWidth', 1.5, 'DisplayName', 'small (legacy)');
plot(t, fM_i_L, 'LineWidth', 1.5, 'DisplayName', 'medium (legacy)');
plot(t, fL_i_L, 'LineWidth', 1.5, 'DisplayName', 'large (legacy)');
plot(t, fS_i_S, '--', 'LineWidth', 1.5, 'DisplayName', 'small (split)');
plot(t, fM_i_S, '--', 'LineWidth', 1.5, 'DisplayName', 'medium (split)');
plot(t, fL_i_S, '--', 'LineWidth', 1.5, 'DisplayName', 'large (split)');
grid on;
xlabel('Time (days)');
ylabel('Export fraction');
title('Export fractions vs time (image diameter)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, sprintf('export_compare_image_%s.png', tag)));

fig2 = figure; hold on;
plot(t, fS_v_L, 'LineWidth', 1.5, 'DisplayName', 'small (legacy)');
plot(t, fM_v_L, 'LineWidth', 1.5, 'DisplayName', 'medium (legacy)');
plot(t, fL_v_L, 'LineWidth', 1.5, 'DisplayName', 'large (legacy)');
plot(t, fS_v_S, '--', 'LineWidth', 1.5, 'DisplayName', 'small (split)');
plot(t, fM_v_S, '--', 'LineWidth', 1.5, 'DisplayName', 'medium (split)');
plot(t, fL_v_S, '--', 'LineWidth', 1.5, 'DisplayName', 'large (split)');
grid on;
xlabel('Time (days)');
ylabel('Export fraction');
title('Export fractions vs time (volume diameter)');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, sprintf('export_compare_volume_%s.png', tag)));

% --------- size spectrum snapshot ---------
snap_t = 100;
[~,idx] = min(abs(t - snap_t));

fig3 = figure; hold on;
loglog(d_cm_v, res_legacy.output_data.masspec_v(idx,:), 'LineWidth', 1.5, 'DisplayName', 'legacy');
loglog(d_cm_v, res_split.output_data.masspec_v(idx,:), '--', 'LineWidth', 1.5, 'DisplayName', 'split');
grid on;
xlabel('Particle diameter [cm]');
ylabel('Mass spectrum [vol/vol/sect/cm]');
title(sprintf('Mass spectrum at t=%.0f d (volume diameter)', t(idx)));
legend('Location','best');
saveas(fig3, fullfile(fig_dir, sprintf('mass_spectrum_compare_%s.png', tag)));

fprintf('\nSaved figures to: %s\n', fig_dir);
