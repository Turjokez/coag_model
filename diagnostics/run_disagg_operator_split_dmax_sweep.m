% run_disagg_operator_split_dmax_sweep.m
% Sweep operator-split Dmax (cm) to see export fraction sensitivity.

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

% Operator-split settings
cfg.disagg_mode = 'operator_split';
cfg.disagg_outer_dt = 1/24;  % days (1 hour)
cfg.disagg_frac_next = 2/3;
cfg.disagg_C = 3.0;
cfg.disagg_gamma = 0.15;

% Dmax sweep (cm)
dmax_list = [2.5, 3, 4, 5, 8];
target_days = [50, 100, 200];

% --------- storage ---------
results = struct();

% --------- run sweep ---------
for i = 1:numel(dmax_list)
    cfg_i = cfg.copy();
    cfg_i.disagg_dmax_cm = dmax_list(i);
    cfg_i.disagg_epsilon = []; % ensure override

    sim = CoagulationSimulation(cfg_i);
    res = sim.run();

    t = res.time(:);
    Fv = res.output_data.fluxspec;
    Fi = res.output_data.fluxspec_i;
    d_cm_v = res.output_data.diam_v(:);
    d_cm_i = res.output_data.diam_i(:);

    fS_v = zeros(numel(t),1); fM_v = fS_v; fL_v = fS_v;
    fS_i = fS_v; fM_i = fS_v; fL_i = fS_v;
    for k = 1:numel(t)
        [fS_v(k), fM_v(k), fL_v(k)] = export_size_classes(d_cm_v, Fv(k,:).', 'volume');
        [fS_i(k), fM_i(k), fL_i(k)] = export_size_classes(d_cm_i, Fi(k,:).', 'image');
    end

    results(i).dmax_cm = dmax_list(i);
    results(i).time = t;
    results(i).fS_v = fS_v;
    results(i).fM_v = fM_v;
    results(i).fL_v = fL_v;
    results(i).fS_i = fS_i;
    results(i).fM_i = fM_i;
    results(i).fL_i = fL_i;

    fprintf('\n=== Dmax = %.2f cm ===\n', dmax_list(i));
    for tt = target_days
        [~,idx] = min(abs(t - tt));
        fprintf('t=%.0f d (vol): S=%.3f M=%.3f L=%.3f | (img): S=%.3f M=%.3f L=%.3f\n', ...
            t(idx), fS_v(idx), fM_v(idx), fL_v(idx), fS_i(idx), fM_i(idx), fL_i(idx));
    end
end

% --------- plots ---------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

colors = lines(numel(dmax_list));

fig1 = figure; hold on;
for i = 1:numel(dmax_list)
    t = results(i).time;
    plot(t, results(i).fS_i, 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('Dmax=%.1f cm', dmax_list(i)));
end
grid on;
xlabel('Time (days)');
ylabel('Small export fraction (image)');
title('Operator-split sweep: small export fraction');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'disagg_split_sweep_small_image.png'));

fig2 = figure; hold on;
for i = 1:numel(dmax_list)
    t = results(i).time;
    plot(t, results(i).fL_v, 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('Dmax=%.1f cm', dmax_list(i)));
end
grid on;
xlabel('Time (days)');
ylabel('Large export fraction (volume)');
title('Operator-split sweep: large export fraction');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'disagg_split_sweep_large_volume.png'));

fprintf('\nSaved sweep figures to: %s\n', fig_dir);
