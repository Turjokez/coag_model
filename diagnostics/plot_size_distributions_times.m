% plot_size_distributions_times.m
% Plot log N vs log D at selected times for each sinking law.
% Saves figures to output/figures.

clear; close all; clc;
setup_paths

laws = {'current','kriest_8','kriest_9','siegel_2025'};

% --------- base configuration ---------
base = SimulationConfig();
base.t_final = 50;
base.delta_t = 1;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;

base.enable_coag = true;
base.enable_sinking = true;
base.enable_disagg = false;
base.enable_linear = true;
base.growth = 0;

base.dz = 65;
base.box_depth = 65;

% times to plot (days)
plot_days = [0 5 10 20 30 40 50];

% output dir
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

for i = 1:numel(laws)
    cfg = copy(base);
    cfg.sinking_law = laws{i};
    cfg.sinking_size = 'volume';
    cfg.sinking_scale = 1.0;

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    t = res.time(:);
    Y = res.concentrations;

    gr = cfg.derive();
    D_mm = gr.getVolumeDiameters() * 10;

    fig = figure('Color','w','Position',[80 80 1100 700]);
    hold on; grid on;

    for k = 1:numel(plot_days)
        [~, idx] = min(abs(t - plot_days(k)));
        yk = Y(idx,:);
        loglog(D_mm, yk, 'DisplayName', sprintf('day %d', plot_days(k)));
    end

    xlabel('Diameter (mm)');
    ylabel('N (state units)');
    title(sprintf('Size distribution (log-log) - %s', laws{i}));
    legend('Location','best');

    out_png = fullfile(fig_dir, sprintf('size_dist_%s.png', laws{i}));
    saveas(fig, out_png);
end

fprintf('Saved size-distribution figures to: %s\n', fig_dir);
