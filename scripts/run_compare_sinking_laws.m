% scripts/run_compare_sinking_laws.m

clear; close all;
setup_paths
clear classes; rehash; rehash toolboxcache

% project root + output folder for figures
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

disp(' ');
disp('Sinking-law comparison (2 laws, 3 figures total)');
disp(['Date: ', datestr(now)]);
disp(' ');

laws = {'kriest_8','kriest_8_capped'};

% Storage
all_total_flux  = [];
all_setvel_mday = [];
all_flux_final  = [];

t    = [];
D_mm = [];

% Fractions (final time)
fS_all = nan(1, numel(laws));
fM_all = nan(1, numel(laws));
fL_all = nan(1, numel(laws));
targets = [50, 100, 200];

% ---------- fixed settings for fair comparison ----------
base = SimulationConfig();
base.t_final = 200;
base.delta_t = 2;
base.n_sections = 35;

base.enable_pp  = true;
base.pp_bin     = 1;
base.pp_source  = 1e-8;

base.alpha = 1.0;

% these switches must be honored by your RHS/LinearProcessBuilder
base.enable_coag   = true;
base.enable_linear = true;
base.enable_sinking = true;
base.enable_disagg = true;
base.c3 = 0.02;
base.r_to_rg = 1.6;

base.growth = 0;

% keep depth consistent
base.dz = 65;
base.box_depth = 2000;   % consistent with main diagnostics

% size-class cutoffs are defined inside export_size_classes.m

for i = 1:numel(laws)

    cfg = copy(base);
    cfg.sinking_law   = laws{i};
    cfg.sinking_size  = 'volume';
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

        gr   = cfg.derive();
        D_mm = gr.getVolumeDiameters() * 10;

    end

    % time series total export flux
    all_total_flux(:,i) = res.output_data.total_flux(:);

    % sinking speed by size (m/day) - should already be in output_data
    all_setvel_mday(:,i) = res.output_data.set_vel(:);

    % flux spectrum final time (volume + image)
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

    % S/M/L fractions using Adrian cutoffs (500/2000 µm)
    d_cm_v = D_mm / 10; % mm -> cm (volume)
    d_cm_i = res.output_data.diam_i(:); % image diameter in cm

    [fS_all(i), fM_all(i), fL_all(i)] = export_size_classes(d_cm_v, flux_final_v, 'volume');
    [fS_i, fM_i, fL_i] = export_size_classes(d_cm_i, flux_final_i, 'image');

    fprintf('%s (vol, final): small=%.3f  med=%.3f  large=%.3f\n', ...
        laws{i}, fS_all(i), fM_all(i), fL_all(i));
    fprintf('%s (img, final): small=%.3f  med=%.3f  large=%.3f\n', ...
        laws{i}, fS_i, fM_i, fL_i);

    for k = 1:numel(targets)
        [~, idx_t] = min(abs(t - targets(k)));
        flux_v_t = Fv(idx_t, :).';
        flux_i_t = Fi(idx_t, :).';
        [sv, mv, lv] = export_size_classes(d_cm_v, flux_v_t, 'volume');
        [si, mi, li] = export_size_classes(d_cm_i, flux_i_t, 'image');
        fprintf('  t=%d d (vol): small=%.3f med=%.3f large=%.3f | (img): %.3f %.3f %.3f\n', ...
            targets(k), sv, mv, lv, si, mi, li);
    end
end

%% Figure 1: total export flux vs time
fig1 = figure; hold on;
for i = 1:numel(laws)
    plot(t, all_total_flux(:,i), 'DisplayName', laws{i});
end
xlabel('Time (days)');
ylabel('Total export flux (model units)');
title('Total export flux vs time');
legend('Location','best'); grid on;
saveas(fig1, fullfile(fig_dir, sprintf('export_flux_c3_%.4g_r%.2f.png', base.c3, base.r_to_rg)));

%% Figure 2: export flux spectrum at final time
fig2 = figure; hold on;
for i = 1:numel(laws)
    loglog(D_mm, all_flux_final(:,i), 'DisplayName', laws{i});
end
xlabel('Diameter (mm)');
ylabel('Export flux by size (model units)');
title('Export flux spectrum at final time');
legend('Location','best'); grid on;
saveas(fig2, fullfile(fig_dir, sprintf('flux_spectrum_c3_%.4g_r%.2f.png', base.c3, base.r_to_rg)));

%% Figure 3: sinking speed vs size
fig3 = figure; hold on;
for i = 1:numel(laws)
    loglog(D_mm, all_setvel_mday(:,i), 'DisplayName', laws{i});
end
xlabel('Diameter (mm)');
ylabel('Sinking speed (m/day)');
title('Sinking speed vs size');
legend('Location','best'); grid on;
saveas(fig3, fullfile(fig_dir, sprintf('sinking_laws_c3_%.4g_r%.2f.png', base.c3, base.r_to_rg)));

%% Print sinking speed near ~1 mm
[~,idx] = min(abs(D_mm - 1.0));
fprintf('\n=== Sinking speed near %.3f mm (bin %d) ===\n', D_mm(idx), idx);
for i = 1:numel(laws)
    fprintf('%s: %.2f m/day\n', laws{i}, all_setvel_mday(idx,i));
end

disp('Done. (3 figures only)');
