% plot_size_distributions_times_image.m
% Plot size distributions at multiple times using image diameter.

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
cfg.enable_disagg = false;
cfg.enable_linear = true;
cfg.growth = 0;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

% --------- run ---------
sim = CoagulationSimulation(cfg);
res = sim.run();

t = res.time(:);
Y = res.concentrations;
d_cm_i = res.output_data.diam_i(:); % image diameter in cm
D_mm_i = d_cm_i * 10;

% Times to plot (nearest)
plot_days = [0, 5, 10, 20, 50, 100, 200];
colors = lines(numel(plot_days));

figure; hold on;
for i = 1:numel(plot_days)
    [~,idx] = min(abs(t - plot_days(i)));
    y = max(Y(idx,:), eps);
    loglog(D_mm_i, y, 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('day %d', round(t(idx))));
end
grid on;
set(gca,'XScale','log','YScale','log');
xlabel('Image diameter (mm)');
ylabel('Concentration (state units)');
title('Size distributions vs time (image diameter)');
legend('Location','best');
% focus on relevant size range
xlim([1e-2 1e2]);
% Adrian cutoffs (500 µm, 2000 µm)
xline(0.5,'--'); xline(2.0,'--');

% Save
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'size_distributions_times_image.png'));

fprintf('Saved figure to: %s\n', fig_dir);
