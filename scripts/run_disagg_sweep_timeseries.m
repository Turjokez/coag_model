% run_disagg_sweep_timeseries.m
% Sweep disaggregation strength and track export fractions vs time.

clear; close all; clc;
setup_paths

c3_vals = [0, 0.02, 0.05, 0.1];

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
base.enable_disagg = true;   % ON for sweep
base.enable_linear = true;
base.growth = 0;

base.sinking_law = 'kriest_8';
base.box_depth = 2000;

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fractions = struct();
t = [];
keys = cell(1, numel(c3_vals));

for i = 1:numel(c3_vals)
    cfg = copy(base);
    cfg.c3 = c3_vals(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    if isempty(t)
        t = res.time(:);
        d_cm_i = res.output_data.diam_i(:);
    end

    F = res.output_data.fluxspec;
    fS = zeros(numel(t),1);
    fM = zeros(numel(t),1);
    fL = zeros(numel(t),1);

    for k = 1:numel(t)
        flux_k = F(k,:).';
        [fS(k), fM(k), fL(k)] = export_size_classes(d_cm_i, flux_k, 'image');
    end

    key_raw = sprintf('c3_%g', c3_vals(i));
    key = matlab.lang.makeValidName(key_raw);
    keys{i} = key;
    fractions.(key) = struct('small', fS, 'med', fM, 'large', fL);

    % Print final time
    fprintf('c3=%g (t=%.0f d): small=%.3f  med=%.3f  large=%.3f\n', ...
        c3_vals(i), t(end), fS(end), fM(end), fL(end));
end

% Plot large fraction vs time
fig1 = figure; hold on;
for i = 1:numel(c3_vals)
    key = keys{i};
    plot(t, fractions.(key).large, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('c3=%g', c3_vals(i)));
end
grid on;
xlabel('Time (days)');
ylabel('Large export fraction (image diam)');
title('Large export fraction vs time (disaggregation sweep)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'disagg_sweep_large_fraction.png'));

% Plot full fractions at final time
fig2 = figure; hold on;
for i = 1:numel(c3_vals)
    key = keys{i};
    plot(t, fractions.(key).small, '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('small c3=%g', c3_vals(i)));
end
grid on;
xlabel('Time (days)');
ylabel('Small export fraction (image diam)');
title('Small export fraction vs time (disaggregation sweep)');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'disagg_sweep_small_fraction.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);
