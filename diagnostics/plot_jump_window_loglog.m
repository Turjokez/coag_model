% plot_jump_window_loglog.m
% Log-log size and flux spectra around the jump window (legacy case).

clear; close all; clc;
setup_paths

cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 0.25;
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = true;
cfg.enable_linear = true;
cfg.growth = 0;
cfg.c3 = 0.02;
cfg.r_to_rg = 1.6;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

sim = CoagulationSimulation(cfg);
res = sim.run();

t = res.time(:);
d_cm_v = res.output_data.diam_v(:);
D_mm_v = d_cm_v * 10;
M_v = max(res.output_data.masspec_v, eps);   % mass spectrum
F_v = max(res.output_data.fluxspec, eps);    % flux spectrum

plot_days = [45 50 55 60 62 64 66 70 80 100];
cmap = turbo(numel(plot_days));

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fig1 = figure('Color','w'); hold on;
for i = 1:numel(plot_days)
    [~, idx] = min(abs(t - plot_days(i)));
    loglog(D_mm_v, M_v(idx,:), 'LineWidth',1.4, ...
        'Color', cmap(i,:), 'DisplayName', sprintf('day %.0f', t(idx)));
end
xline(0.5, 'k--', '500 um');
xline(2.0, 'k--', '2000 um');
grid on;
xlabel('Diameter (mm)');
ylabel('Mass spectrum');
title('Legacy mass spectrum (log-log, jump window)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'legacy_jump_window_mass_loglog.png'));

fig2 = figure('Color','w'); hold on;
for i = 1:numel(plot_days)
    [~, idx] = min(abs(t - plot_days(i)));
    loglog(D_mm_v, F_v(idx,:), 'LineWidth',1.4, ...
        'Color', cmap(i,:), 'DisplayName', sprintf('day %.0f', t(idx)));
end
xline(0.5, 'k--', '500 um');
xline(2.0, 'k--', '2000 um');
grid on;
xlabel('Diameter (mm)');
ylabel('Flux spectrum');
title('Legacy flux spectrum (log-log, jump window)');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'legacy_jump_window_flux_loglog.png'));

fprintf('\nSaved jump-window log-log figures to: %s\n', fig_dir);

