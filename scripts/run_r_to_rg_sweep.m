% run_r_to_rg_sweep.m
% Sweep r_to_rg (image scaling) and report image-based export fractions.

clear; close all; clc;
setup_paths

rvals = [1.0, 1.2, 1.4, 1.6, 1.8];

% --------- base configuration ---------
base = SimulationConfig();
base.t_final = 200;
base.delta_t = 2;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;

base.enable_coag = true;
base.enable_sinking = true;
base.enable_disagg = true;
base.c3 = 0.02;
base.enable_linear = true;
base.growth = 0;

base.sinking_law = 'kriest_8';
base.sinking_size = 'volume';
base.box_depth = 2000;

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fS = nan(numel(rvals), 3);
fM = nan(numel(rvals), 3);
fL = nan(numel(rvals), 3);
targets = [50, 100, 200];

for i = 1:numel(rvals)
    cfg = copy(base);
    cfg.r_to_rg = rvals(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    t = res.time(:);
    d_cm_i = res.output_data.diam_i(:);
    Fi = res.output_data.fluxspec_i;
    if size(Fi,1) ~= numel(t)
        Fi = Fi';
    end

    for k = 1:numel(targets)
        [~, idx] = min(abs(t - targets(k)));
        flux_i = Fi(idx, :).';
        [fS(i, k), fM(i, k), fL(i, k)] = export_size_classes(d_cm_i, flux_i, 'image');
    end

    fprintf('r_to_rg=%.2f: t=50 (S/M/L)=%.3f/%.3f/%.3f | t=100 %.3f/%.3f/%.3f | t=200 %.3f/%.3f/%.3f\n', ...
        rvals(i), fS(i,1), fM(i,1), fL(i,1), fS(i,2), fM(i,2), fL(i,2), fS(i,3), fM(i,3), fL(i,3));
end

% Plot fractions vs r_to_rg for each target time
for k = 1:numel(targets)
    fig1 = figure; hold on;
    plot(rvals, fS(:,k), '-o', 'DisplayName', 'small');
    plot(rvals, fM(:,k), '-o', 'DisplayName', 'medium');
    plot(rvals, fL(:,k), '-o', 'DisplayName', 'large');
    grid on;
    xlabel('r\_to\_rg');
    ylabel('Export fraction (image)');
    title(sprintf('Export fractions vs r\\_to\\_rg (t=%d d, c3=%.4g)', targets(k), base.c3));
    legend('Location','best');
    saveas(fig1, fullfile(fig_dir, sprintf('r_to_rg_sweep_export_fractions_t%d_c3_%.4g.png', targets(k), base.c3)));
end

fprintf('\nSaved figure to: %s\n', fig_dir);
