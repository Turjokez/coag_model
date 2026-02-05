% run_kernel_scale_sweep.m
% Sweep alpha over orders of magnitude with NO sinking and zero init.
% Purpose: check if coagulation alone can build large sizes.

clear; close all; clc;
setup_paths

alphas = [1, 10, 1e2, 1e3, 1e4];

% --------- base configuration ---------
base = SimulationConfig();
base.t_final = 200;
base.delta_t = 2;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;

base.enable_coag = true;
base.enable_sinking = false;  % OFF
base.enable_disagg = false;
base.enable_linear = true;
base.growth = 0;

base.sinking_law = 'kriest_8';
base.box_depth = 2000;

v0 = zeros(base.n_sections,1);

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

D_mm = [];
spec_all = [];
frac_small = nan(1, numel(alphas));
frac_med   = nan(1, numel(alphas));
frac_large = nan(1, numel(alphas));

for i = 1:numel(alphas)
    cfg = copy(base);
    cfg.alpha = alphas(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run('v0', v0);

    if isempty(D_mm)
        gr = cfg.derive();
        D_mm = gr.getVolumeDiameters() * 10;
    end

    v_end = res.concentrations(end,:)';
    v_end = max(v_end, 0);
    spec_all(:,i) = v_end; %#ok<AGROW>

    % size-class fractions based on concentration (not export)
    d_cm = D_mm / 10;
    idx_small = d_cm < 0.05;              % < 500 µm
    idx_med   = d_cm >= 0.05 & d_cm < 0.2; % 0.5–2 mm
    idx_large = d_cm >= 0.2;              % > 2 mm

    total = sum(v_end);
    if total > 0
        frac_small(i) = sum(v_end(idx_small)) / total;
        frac_med(i)   = sum(v_end(idx_med)) / total;
        frac_large(i) = sum(v_end(idx_large)) / total;
    else
        frac_small(i) = 0; frac_med(i) = 0; frac_large(i) = 0;
    end

    fprintf('alpha=%g: small=%.3f  med=%.3f  large=%.3f (concentration fractions)\n', ...
        alphas(i), frac_small(i), frac_med(i), frac_large(i));
end

% Plot size distributions
fig1 = figure; hold on;
for i = 1:numel(alphas)
    loglog(D_mm, max(spec_all(:,i), eps), 'DisplayName', sprintf('alpha=%g', alphas(i)));
end
grid on;
xlabel('Diameter (mm)');
ylabel('Concentration (state units)');
title('Size distribution (day 200, zero init, no sinking)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'kernel_scale_size_distributions.png'));

% Plot fractions vs alpha
fig2 = figure; hold on;
semilogx(alphas, frac_small, '-o', 'DisplayName', 'small');
semilogx(alphas, frac_med, '-o', 'DisplayName', 'medium');
semilogx(alphas, frac_large, '-o', 'DisplayName', 'large');
grid on;
xlabel('alpha (stickiness)');
ylabel('Concentration fraction');
title('Size-class fractions vs alpha (no sinking)');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'kernel_scale_fractions.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);
