% plot_rate_terms_by_size.m
% Plot rate terms by size at selected times (includes disagg).

clear; close all; clc;
setup_paths

cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 2;
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = true;
cfg.c3 = 0.02;
cfg.r_to_rg = 1.6;
cfg.enable_linear = true;
cfg.growth = 0;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

sim = CoagulationSimulation(cfg);
res = sim.run();

plot_days = [50 100 200];

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

tag = sprintf('c3_%.4g_r%.2f', cfg.c3, cfg.r_to_rg);
MassBalanceBiovolume.plotRateTermsBySize(sim, res, plot_days, fig_dir, tag);

fprintf('Saved rate-term plots to: %s\n', fig_dir);
