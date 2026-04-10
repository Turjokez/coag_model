% build_apr08_1d_budget_split
% Short note:
% 1. rerun a few trusted 1-D cases
% 2. plot volume in column, volume out bottom, tracked total
% 3. save one simple panel for the report

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
work_root = fileparts(repo_root);
test_root = fullfile(work_root, '1d-model-testing');
addpath(genpath(fullfile(test_root, 'src')));

fig_dir = fullfile(repo_root, 'docs', 'figures');

sim_diff = build_step3_case();
sim_coag = build_step4_case();
sim_frag = build_step5_case();
sim_sink = build_step7_case();

fig = figure('Color', 'w', 'Position', [120 120 1500 1180]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_budget_panel(nexttile(tl, 1), sim_diff, 'Step 3: diffusion');
plot_budget_panel(nexttile(tl, 2), sim_coag, 'Step 4: coagulation');
plot_budget_panel(nexttile(tl, 3), sim_frag, 'Step 5: fragmentation');
plot_budget_panel(nexttile(tl, 4), sim_sink, 'Step 7: sinking coupling');

save_figure(fig, fullfile(fig_dir, 'apr08_1d_panel_budget_split.png'));
close(fig);

disp('Saved budget split figure:');
disp(fullfile(fig_dir, 'apr08_1d_panel_budget_split.png'));

function sim = build_step3_case()
law_name = 'kriest_8';
size_um = [100; 500; 1000; 3000];
pulse_amp = [1.0; 0.6; 0.3; 0.15];
size_cm = size_um .* 1e-4;
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(speed_m_s);
cfg.t_max_s = 1.20 .* max(cfg.z_max_m ./ speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;

sim = solve_advection_diffusion(cfg);
end

function sim = build_step4_case()
law_name = 'kriest_8';
size_um = round(logspace(log10(200), log10(3000), 8))';
size_cm = size_um .* 1e-4;
pulse_amp = powerlaw_concentration(size_cm, 5e-3, -2.5);
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(speed_m_s);
cfg.t_max_s = 1.10 .* max(cfg.z_max_m ./ speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;
cfg.kernel_mode = 'shear_only';
cfg.epsilon_mks = 1e-6;
cfg.coag_scale = 100.0;
cfg.coag_substeps = 4;
cfg.scale_shear = 1.0;
cfg.scale_diff_sed = 0.0;

sim = solve_with_coagulation(cfg);
end

function sim = build_step5_case()
law_name = 'kriest_8';
size_um = round(logspace(log10(200), log10(3000), 8))';
size_cm = size_um .* 1e-4;
pulse_amp = powerlaw_concentration(size_cm, 5e-3, -2.5);
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(speed_m_s);
cfg.t_max_s = 1.10 .* max(cfg.z_max_m ./ speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;
cfg.kernel_mode = 'shear_only';
cfg.epsilon_mks = 1e-6;
cfg.coag_scale = 100.0;
cfg.coag_substeps = 4;
cfg.scale_shear = 1.0;
cfg.scale_diff_sed = 0.0;
cfg.frag_substeps = 4;
cfg.c3 = 0.005;
cfg.c4 = 1.45;

sim = solve_with_fragmentation(cfg);
end

function sim = build_step7_case()
law_name = 'kriest_8';
size_um = round(logspace(log10(200), log10(3000), 8))';
size_cm = size_um .* 1e-4;
pulse_amp = powerlaw_concentration(size_cm, 5e-3, -2.5);
base_speed_cm_s = sinking_speed_named(size_cm, law_name);
base_speed_m_s = base_speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(base_speed_m_s);
cfg.t_max_s = 1.10 .* max(cfg.z_max_m ./ base_speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = base_speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;
cfg.kernel_mode = 'shear_only';
cfg.epsilon_mks = 1e-6;
cfg.coag_scale = 100.0;
cfg.coag_substeps = 4;
cfg.scale_shear = 1.0;
cfg.scale_diff_sed = 0.0;
cfg.frag_substeps = 4;
cfg.c3 = 0.005;
cfg.c4 = 1.45;

grid = make_depth_grid(cfg.z_max_m, cfg.dz_m);
prof = local_profiles(grid.z_m);
cfg.kz_profile_m2_s = prof.kz_m2_s;
cfg.temp_profile_c = prof.temp_c;
cfg.sal_profile_psu = prof.sal_psu;
cfg.rho_profile_kg_m3 = prof.rho_kg_m3;

sink_prof = build_sinking_speed_profile(base_speed_m_s, prof.temp_c, prof.rho_kg_m3);
cfg.speed_profile_m_s = sink_prof.speed_profile_m_s;

sim = solve_with_fragmentation(cfg);
end

function plot_budget_panel(ax, sim, ttl)
t_day = sim.t_s ./ 86400.0;
v0 = max(sim.tracked_volume_total(1), realmin);
in_col = 100.0 .* sim.column_volume_total ./ v0;
out_bot = 100.0 .* sim.export_volume_total ./ v0;
track = 100.0 .* sim.tracked_volume_total ./ v0;

hold(ax, 'on');
plot(ax, t_day, in_col, 'k', 'LineWidth', 1.5, 'DisplayName', 'volume in column');
plot(ax, t_day, out_bot, 'r', 'LineWidth', 1.5, 'DisplayName', 'volume out bottom');
plot(ax, t_day, track, 'b', 'LineWidth', 1.5, 'DisplayName', 'tracked total');
xlabel(ax, 'Time (day)');
ylabel(ax, 'Volume (% of start)');
title(ax, ttl);
ylim(ax, [0 105]);

if strcmp(ttl, 'Step 3: diffusion')
    legend(ax, 'Location', 'best', 'Box', 'off');
end

ax.LineWidth = 1.0;
ax.FontSize = 11;
end

function prof = local_profiles(z_m)
prof = struct();
prof.temp_c = 4.0 + 14.0 .* exp(-z_m ./ 150.0);
prof.sal_psu = 34.2 + 0.8 .* (1.0 - exp(-z_m ./ 300.0));
prof.rho_kg_m3 = 1024.5 + 2.2 .* (1.0 - exp(-z_m ./ 250.0));
prof.kz_m2_s = 1e-6 + (1.5e-3 - 1e-6) .* exp(-z_m ./ 120.0);
end
