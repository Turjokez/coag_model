% run_kernel_component_sweep.m
% Compare kernel components by turning each on/off.

clear; close all; clc;
setup_paths

cases = {
    'all',     1, 1, 1;
    'brown',   1, 0, 0;
    'shear',   0, 1, 0;
    'ds',      0, 0, 1;
    'no_ds',   1, 1, 0;
};

% --------- base configuration ---------
base = SimulationConfig();
base.t_final = 200;
base.delta_t = 2;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;

base.enable_coag = true;
base.enable_sinking = false;  % OFF for pure coag test
base.enable_disagg = false;
base.enable_linear = true;
base.growth = 0;

base.sinking_law = 'kriest_8';
base.box_depth = 2000;
base.alpha = 1.0;

v0 = zeros(base.n_sections,1);

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

D_mm = [];
spec_all = [];
labels = {};

for i = 1:size(cases,1)
    name = cases{i,1};
    sb = cases{i,2};
    ss = cases{i,3};
    sd = cases{i,4};

    cfg = copy(base);
    cfg.scale_brown = sb;
    cfg.scale_shear = ss;
    cfg.scale_ds    = sd;

    sim = CoagulationSimulation(cfg);
    res = sim.run('v0', v0);

    if isempty(D_mm)
        gr = cfg.derive();
        D_mm = gr.getVolumeDiameters() * 10;
    end

    v_end = max(res.concentrations(end,:)', 0);
    spec_all(:,i) = v_end; %#ok<AGROW>
    labels{i} = name; %#ok<AGROW>

    % size-class fractions based on concentration
    d_cm = D_mm / 10;
    idx_small = d_cm < 0.05;
    idx_med   = d_cm >= 0.05 & d_cm < 0.2;
    idx_large = d_cm >= 0.2;
    total = sum(v_end);
    if total > 0
        fS = sum(v_end(idx_small)) / total;
        fM = sum(v_end(idx_med)) / total;
        fL = sum(v_end(idx_large)) / total;
    else
        fS = 0; fM = 0; fL = 0;
    end
    fprintf('%s: small=%.3f  med=%.3f  large=%.3f (concentration fractions)\n', ...
        name, fS, fM, fL);
end

% Plot size distributions
fig1 = figure; hold on;
for i = 1:numel(labels)
    loglog(D_mm, max(spec_all(:,i), eps), 'DisplayName', labels{i});
end
grid on;
xlabel('Diameter (mm)');
ylabel('Concentration (state units)');
title('Size distribution (day 200, zero init, no sinking)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'kernel_component_size_distributions.png'));

fprintf('\nSaved figure to: %s\n', fig_dir);
