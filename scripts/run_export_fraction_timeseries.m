% run_export_fraction_timeseries.m
% Track export size-class fractions over time (volume vs image diameter).

clear; close all; clc;
setup_paths

% --------- configuration ---------
cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 2;
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = true;
cfg.c3 = 0.02;
cfg.r_to_rg = 1.6;
cfg.enable_linear = true;
cfg.growth = 0;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

% --------- run ---------
sim = CoagulationSimulation(cfg);
res = sim.run();

t = res.time(:);
Fv = res.output_data.fluxspec;   % time x bins (volume)
Fi = res.output_data.fluxspec_i; % time x bins (image)

d_cm_v = res.output_data.diam_v(:);
d_cm_i = res.output_data.diam_i(:);

fS_v = zeros(numel(t),1);
fM_v = zeros(numel(t),1);
fL_v = zeros(numel(t),1);

fS_i = zeros(numel(t),1);
fM_i = zeros(numel(t),1);
fL_i = zeros(numel(t),1);

for k = 1:numel(t)
    flux_k_v = Fv(k,:).';
    flux_k_i = Fi(k,:).';
    [fS_v(k), fM_v(k), fL_v(k)] = export_size_classes(d_cm_v, flux_k_v, 'volume');
    [fS_i(k), fM_i(k), fL_i(k)] = export_size_classes(d_cm_i, flux_k_i, 'image');
end

% Print fractions at specific times (nearest)
targets = [50, 100, 200];
fprintf('\n=== Export fractions (nearest times) ===\n');
for tt = targets
    [~,idx] = min(abs(t - tt));
    fprintf('t=%.0f d (vol): small=%.3f med=%.3f large=%.3f\n', ...
        t(idx), fS_v(idx), fM_v(idx), fL_v(idx));
    fprintf('t=%.0f d (img): small=%.3f med=%.3f large=%.3f\n', ...
        t(idx), fS_i(idx), fM_i(idx), fL_i(idx));
end

% --------- plots ---------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
tag = sprintf('c3_%.4g', cfg.c3);

fig1 = figure; hold on;
plot(t, fS_i, 'LineWidth', 1.5, 'DisplayName', 'small (image)');
plot(t, fM_i, 'LineWidth', 1.5, 'DisplayName', 'medium (image)');
plot(t, fL_i, 'LineWidth', 1.5, 'DisplayName', 'large (image)');
grid on;
xlabel('Time (days)');
ylabel('Export fraction');
title(sprintf('Export fractions vs time (image diameter), c3=%.4g', cfg.c3));
legend('Location','best');
saveas(fig1, fullfile(fig_dir, sprintf('export_fractions_time_image_%s.png', tag)));

fig2 = figure; hold on;
plot(t, fS_v, 'LineWidth', 1.5, 'DisplayName', 'small (volume)');
plot(t, fM_v, 'LineWidth', 1.5, 'DisplayName', 'medium (volume)');
plot(t, fL_v, 'LineWidth', 1.5, 'DisplayName', 'large (volume)');
grid on;
xlabel('Time (days)');
ylabel('Export fraction');
title(sprintf('Export fractions vs time (volume diameter), c3=%.4g', cfg.c3));
legend('Location','best');
saveas(fig2, fullfile(fig_dir, sprintf('export_fractions_time_volume_%s.png', tag)));

fprintf('\nSaved figures to: %s\n', fig_dir);
