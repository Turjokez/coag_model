% run_pp_source_sweep.m
% Sweep PP source strength and compare export fractions + spectra.

clear; close all; clc;
setup_paths

pp_sources = [1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5];

% --------- base configuration ---------
base = SimulationConfig();
base.t_final = 50;
base.delta_t = 1;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;

base.enable_coag = true;
base.enable_sinking = true;
base.enable_disagg = false;
base.enable_linear = true;
base.growth = 0;

base.sinking_law = 'kriest_8';
base.sinking_size = 'volume';
base.sinking_scale = 1.0;
base.box_depth = 2000;
base.alpha = 1.0;

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fS_all = nan(1, numel(pp_sources));
fM_all = nan(1, numel(pp_sources));
fL_all = nan(1, numel(pp_sources));

D_mm = [];
flux_final_all = [];

for i = 1:numel(pp_sources)
    cfg = copy(base);
    cfg.pp_source = pp_sources(i);

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
    fprintf('pp_source=%g: small=%.3f  med=%.3f  large=%.3f\n', ...
        pp_sources(i), fS_all(i), fM_all(i), fL_all(i));
end

% Figure 1: export spectrum at day 50
fig1 = figure; hold on;
for i = 1:numel(pp_sources)
    loglog(D_mm, flux_final_all(:,i), 'DisplayName', sprintf('PP=%g', pp_sources(i)));
end
grid on;
xlabel('Diameter (mm)');
ylabel('Export flux by size (model units)');
title('Export flux spectrum at day 50 (PP sweep)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'pp_sweep_flux_spectrum.png'));

% Figure 2: S/M/L fractions vs PP
fig2 = figure; hold on;
semilogx(pp_sources, fS_all, '-o', 'DisplayName', 'small');
semilogx(pp_sources, fM_all, '-o', 'DisplayName', 'medium');
semilogx(pp_sources, fL_all, '-o', 'DisplayName', 'large');
grid on;
xlabel('PP source (state units / day)');
ylabel('Export fraction');
title('Export fractions vs PP source');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'pp_sweep_export_fractions.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);
