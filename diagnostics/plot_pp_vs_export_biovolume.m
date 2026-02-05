% plot_pp_vs_export_biovolume.m
% Compare PP input vs sinking export in consistent (state/biovolume) units.

clear; close all; clc;
setup_paths

% --------- configuration (edit if needed) ---------
cfg = SimulationConfig();
cfg.t_final = 50;
cfg.delta_t = 1;
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = false;
cfg.enable_linear = true;
cfg.growth = 0;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

% --------- run ---------
sim = CoagulationSimulation(cfg);
res = sim.run();

t = res.time(:);
Y = res.concentrations;

% PP input rate from RHS (state units / day)
pp_rate = zeros(numel(t),1);
for i = 1:numel(t)
    [~,~,~,~,term5] = sim.rhs.rateTerms(Y(i,:)');
    pp_rate(i) = sum(term5);
end

% Sinking export loss (state units / day)
export_rate = res.output_data.sinkLossTotal(:);

% --------- plot ---------
figure; hold on;
plot(t, pp_rate, 'LineWidth', 1.5, 'DisplayName', 'PP input');
plot(t, export_rate, 'LineWidth', 1.5, 'DisplayName', 'Sinking export');
grid on;
xlabel('Time (days)');
ylabel('Rate (state units / day)');
title('PP input vs sinking export (biovolume-consistent)');
legend('Location','best');

fprintf('PP input (constant): %.3e state units/day\n', pp_rate(1));
fprintf('Export rate (final): %.3e state units/day\n', export_rate(end));
