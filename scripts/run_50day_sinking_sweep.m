% run_50day_sinking_sweep.m
% 50-day sweep of 4 sinking laws with consistent setup.
% Saves figures to output/figures and prints S/M/L export fractions.

clear; close all; clc;
setup_paths

laws = {'current','kriest_8','kriest_8_capped','kriest_8_flat','kriest_9','siegel_2025'};

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
base.growth = 0; % keep PP separate for budget clarity

% keep depth consistent
base.dz = 65;
base.box_depth = 65;

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% storage
all_total_flux = [];
all_setvel_mday = [];
all_flux_final = [];

fS_all = nan(1, numel(laws));
fM_all = nan(1, numel(laws));
fL_all = nan(1, numel(laws));

t = [];
D_mm = [];

for i = 1:numel(laws)
    cfg = copy(base);
    cfg.sinking_law = laws{i};
    cfg.sinking_size = 'volume';
    cfg.sinking_scale = 1.0;
    if strcmp(cfg.sinking_law, 'kriest_8_capped')
        cfg.sinking_w_max_mday = 70; % cap per Adrian-style diagnostic
    end
    if strcmp(cfg.sinking_law, 'kriest_8_flat')
        cfg.sinking_d_flat_cm = 0.1; % 1 mm flat threshold
    end

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    if isempty(t)
        t = res.time(:);
        gr = cfg.derive();
        D_mm = gr.getVolumeDiameters() * 10;
    end

    % time series total export flux
    all_total_flux(:,i) = res.output_data.total_flux(:);

    % sinking speed by size (m/day)
    all_setvel_mday(:,i) = res.output_data.set_vel(:);

    % flux spectrum at final time (volume + image)
    Fv = res.output_data.fluxspec;
    Fi = res.output_data.fluxspec_i;
    if size(Fv,1) == numel(t)
        flux_final_v = Fv(end,:).';
        flux_final_i = Fi(end,:).';
    else
        flux_final_v = Fv(:,end);
        flux_final_i = Fi(:,end);
    end
    all_flux_final(:,i) = flux_final_v;

    % S/M/L fractions using export_size_classes (volume + image)
    d_cm_v = D_mm / 10; % mm -> cm (volume)
    d_cm_i = res.output_data.diam_i(:); % image diameter in cm

    [fS_all(i), fM_all(i), fL_all(i)] = export_size_classes(d_cm_v, flux_final_v, 'volume');
    [fS_i, fM_i, fL_i] = export_size_classes(d_cm_i, flux_final_i, 'image');

    fprintf('%s (vol): small=%.3f  med=%.3f  large=%.3f\n', ...
        laws{i}, fS_all(i), fM_all(i), fL_all(i));
    fprintf('%s (img): small=%.3f  med=%.3f  large=%.3f\n', ...
        laws{i}, fS_i, fM_i, fL_i);
end

%% Figure 1: total export flux vs time
fig1 = figure; hold on;
for i = 1:numel(laws)
    plot(t, all_total_flux(:,i), 'DisplayName', laws{i});
end
xlabel('Time (days)');
ylabel('Total export flux (model units)');
title('Total export flux vs time (50-day sweep)');
legend('Location','best'); grid on;
saveas(fig1, fullfile(fig_dir, 'sweep50_export_flux.png'));

%% Figure 2: export flux spectrum at final time
fig2 = figure; hold on;
for i = 1:numel(laws)
    loglog(D_mm, all_flux_final(:,i), 'DisplayName', laws{i});
end
xlabel('Diameter (mm)');
ylabel('Export flux by size (model units)');
title('Export flux spectrum at day 50');
legend('Location','best'); grid on;
saveas(fig2, fullfile(fig_dir, 'sweep50_flux_spectrum.png'));

%% Figure 3: sinking speed vs size
fig3 = figure; hold on;
for i = 1:numel(laws)
    loglog(D_mm, all_setvel_mday(:,i), 'DisplayName', laws{i});
end
xlabel('Diameter (mm)');
ylabel('Sinking speed (m/day)');
title('Sinking speed vs size');
legend('Location','best'); grid on;
saveas(fig3, fullfile(fig_dir, 'sweep50_sinking_laws.png'));

%% Print sinking speed near ~1 mm
[~,idx] = min(abs(D_mm - 1.0));
fprintf('\n=== Sinking speed near %.3f mm (bin %d) ===\n', D_mm(idx), idx);
for i = 1:numel(laws)
    fprintf('%s: %.2f m/day\n', laws{i}, all_setvel_mday(idx,i));
end

fprintf('\nSaved figures to: %s\n', fig_dir);
