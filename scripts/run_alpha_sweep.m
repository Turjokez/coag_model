% run_alpha_sweep.m
% Sweep coagulation stickiness alpha and compare export fractions.

clear; close all; clc;
setup_paths

alphas = [1, 2, 5, 10, 20];

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

base.sinking_law = 'kriest_8';
base.sinking_size = 'volume';
base.sinking_scale = 1.0;
base.box_depth = 2000;

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fS_all = nan(1, numel(alphas));
fM_all = nan(1, numel(alphas));
fL_all = nan(1, numel(alphas));

D_mm = [];
flux_final_all = [];

for i = 1:numel(alphas)
    cfg = copy(base);
    cfg.alpha = alphas(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    if isempty(D_mm)
        gr = cfg.derive();
        D_mm = gr.getVolumeDiameters() * 10;
    end

    F = res.output_data.fluxspec;
    if size(F,1) == numel(res.time)
        flux_final = F(end,:).';
    else
        flux_final = F(:,end);
    end
    flux_final_all(:,i) = flux_final; %#ok<AGROW>

    [fS_all(i), fM_all(i), fL_all(i)] = export_size_classes(D_mm/10, flux_final);
    fprintf('alpha=%g: small=%.3f  med=%.3f  large=%.3f\n', ...
        alphas(i), fS_all(i), fM_all(i), fL_all(i));
end

% Plot export spectrum at day 50 for each alpha
fig1 = figure; hold on;
for i = 1:numel(alphas)
    loglog(D_mm, flux_final_all(:,i), 'DisplayName', sprintf('alpha=%g', alphas(i)));
end
grid on;
xlabel('Diameter (mm)');
ylabel('Export flux by size (model units)');
title('Export flux spectrum at day 50 (alpha sweep)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'alpha_sweep_flux_spectrum.png'));

fprintf('\nSaved figure to: %s\n', fig_dir);
