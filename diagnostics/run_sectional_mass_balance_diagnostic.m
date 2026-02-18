% run_sectional_mass_balance_diagnostic.m
% Detailed sectional in/out diagnostics for legacy case.

clear; close all; clc;
setup_paths

cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 0.25;
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = true;
cfg.enable_linear = true;
cfg.growth = 0;
cfg.c3 = 0.02;
cfg.r_to_rg = 1.6;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

sim = CoagulationSimulation(cfg);
res = sim.run();
mb = SectionalMassBalanceDetailed(sim, res);

t = mb.t;
d_cm = res.output_data.diam_v(:);
n_sec = numel(d_cm);

[~, i05] = min(abs(d_cm - 0.05));
[~, i20] = min(abs(d_cm - 0.20));
sel = unique(max(1, min(n_sec, [1, i05, i20, n_sec])));

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% --------- closure summary ---------
max_abs_fd = max(abs(mb.closure_fd), [], 2);
max_abs_rhs = max(abs(mb.closure_rhs), [], 2);

fig1 = figure('Color','w');
subplot(2,1,1);
plot(t, max_abs_fd, 'LineWidth',1.5);
grid on; xlabel('Time (days)'); ylabel('max |closure fd|');
title('Sectional closure: finite-difference check');

subplot(2,1,2);
plot(t, max_abs_rhs, 'LineWidth',1.5);
grid on; xlabel('Time (days)'); ylabel('max |closure rhs|');
title('Sectional closure: RHS check');
saveas(fig1, fullfile(fig_dir, 'sectional_balance_closure_summary.png'));

fig2 = figure('Color','w');
imagesc(t, 1:n_sec, mb.closure_fd');
axis xy; colorbar;
xlabel('Time (days)');
ylabel('Section index');
title('Sectional closure residual (dY/dt - net in/out)');
saveas(fig2, fullfile(fig_dir, 'sectional_balance_closure_heatmap.png'));

% --------- total in/out for selected bins ---------
fig3 = figure('Color','w');
for k = 1:numel(sel)
    ib = sel(k);
    subplot(numel(sel),1,k); hold on;
    plot(t, mb.total_in(:,ib), 'LineWidth',1.2, 'DisplayName','total in');
    plot(t, mb.total_out(:,ib), 'LineWidth',1.2, 'DisplayName','total out');
    plot(t, mb.net_model(:,ib), 'LineWidth',1.2, 'DisplayName','net model');
    plot(t, mb.dYdt_fd(:,ib), '--', 'LineWidth',1.1, 'DisplayName','dY/dt fd');
    grid on;
    ylabel(sprintf('bin %d', ib));
    title(sprintf('D=%.3f cm', d_cm(ib)));
    if k == 1, legend('Location','best'); end
end
xlabel('Time (days)');
saveas(fig3, fullfile(fig_dir, 'sectional_balance_selected_bins_total.png'));

% --------- process-wise in/out for boundary bins ---------
fig4 = figure('Color','w');
bins_proc = unique([i05 i20]);
for k = 1:numel(bins_proc)
    ib = bins_proc(k);
    subplot(numel(bins_proc),1,k); hold on;
    plot(t, mb.coag_in(:,ib), 'LineWidth',1.2, 'DisplayName','coag in');
    plot(t, mb.coag_out(:,ib), 'LineWidth',1.2, 'DisplayName','coag out');
    plot(t, mb.growth_in(:,ib), 'LineWidth',1.2, 'DisplayName','growth in');
    plot(t, mb.growth_out(:,ib), 'LineWidth',1.2, 'DisplayName','growth out');
    plot(t, mb.sink_out(:,ib), 'LineWidth',1.2, 'DisplayName','sink out');
    plot(t, mb.disagg_in(:,ib), 'LineWidth',1.2, 'DisplayName','disagg in');
    plot(t, mb.disagg_out(:,ib), 'LineWidth',1.2, 'DisplayName','disagg out');
    plot(t, mb.pp_in(:,ib), 'LineWidth',1.2, 'DisplayName','pp in');
    grid on;
    ylabel(sprintf('bin %d', ib));
    title(sprintf('Process terms at D=%.3f cm', d_cm(ib)));
    if k == 1, legend('Location','best'); end
end
xlabel('Time (days)');
saveas(fig4, fullfile(fig_dir, 'sectional_balance_boundary_bins_process_terms.png'));

fprintf('\nSaved sectional mass-balance figures to: %s\n', fig_dir);

