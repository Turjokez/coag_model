% run_legacy_export_jump_diagnostic.m
% Debug the sharp transition in legacy export fractions.

clear; close all; clc;
setup_paths

% --------- config (legacy) ---------
cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 0.25;    % finer output near jump
cfg.n_sections = 35;

cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;

cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_disagg = true;   % legacy disagg in RHS
cfg.enable_linear = true;
cfg.growth = 0;
cfg.c3 = 0.02;
cfg.r_to_rg = 1.6;

cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;

sim = CoagulationSimulation(cfg);
res = sim.run();

t = res.time(:);
Fv = max(res.output_data.fluxspec, 0);
Fi = max(res.output_data.fluxspec_i, 0);
d_cm_v = res.output_data.diam_v(:);
d_cm_i = res.output_data.diam_i(:);

% --------- class fractions + absolute export ---------
fS_v = zeros(numel(t),1); fM_v = fS_v; fL_v = fS_v;
fS_i = fS_v; fM_i = fS_v; fL_i = fS_v;
S_v = fS_v; M_v = fS_v; L_v = fS_v; T_v = fS_v;

for k = 1:numel(t)
    [fS_v(k), fM_v(k), fL_v(k), tot_v] = export_size_classes(d_cm_v, Fv(k,:).', 'volume');
    [fS_i(k), fM_i(k), fL_i(k)] = export_size_classes(d_cm_i, Fi(k,:).', 'image');
    S_v(k) = tot_v.small;
    M_v(k) = tot_v.medium;
    L_v(k) = tot_v.large;
    T_v(k) = tot_v.total;
end

% detect strongest change point in legacy volume fractions
df = abs(diff(fS_v)) + abs(diff(fM_v)) + abs(diff(fL_v));
[~, ix_jump] = max(df);
t_jump = t(ix_jump+1);

fprintf('\n=== Legacy jump diagnostic ===\n');
fprintf('Largest fraction change near t = %.2f days\n', t_jump);
fprintf('Before (t=%.2f): S=%.3f M=%.3f L=%.3f, total=%.3e\n', ...
    t(ix_jump), fS_v(ix_jump), fM_v(ix_jump), fL_v(ix_jump), T_v(ix_jump));
fprintf('After  (t=%.2f): S=%.3f M=%.3f L=%.3f, total=%.3e\n', ...
    t(ix_jump+1), fS_v(ix_jump+1), fM_v(ix_jump+1), fL_v(ix_jump+1), T_v(ix_jump+1));

% --------- boundary-bin contribution check ---------
[~, i05] = min(abs(d_cm_v - 0.05));
[~, i20] = min(abs(d_cm_v - 0.20));

bins05 = max(1, i05-2):min(numel(d_cm_v), i05+2);
bins20 = max(1, i20-2):min(numel(d_cm_v), i20+2);

win = (t >= max(min(t), t_jump-15)) & (t <= min(max(t), t_jump+15));
t_win = t(win);

% --------- save plots ---------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

fig1 = figure('Color','w');
subplot(2,1,1); hold on;
plot(t, fS_v, 'LineWidth',1.4, 'DisplayName','small (vol)');
plot(t, fM_v, 'LineWidth',1.4, 'DisplayName','medium (vol)');
plot(t, fL_v, 'LineWidth',1.4, 'DisplayName','large (vol)');
plot(t, fS_i, '--', 'LineWidth',1.2, 'DisplayName','small (img)');
plot(t, fM_i, '--', 'LineWidth',1.2, 'DisplayName','medium (img)');
plot(t, fL_i, '--', 'LineWidth',1.2, 'DisplayName','large (img)');
xline(t_jump, 'k--', 'LineWidth',1.2, 'DisplayName','jump');
grid on;
xlabel('Time (days)'); ylabel('Export fraction');
title('Legacy export fractions');
legend('Location','best');

subplot(2,1,2); hold on;
plot(t, S_v, 'LineWidth',1.4, 'DisplayName','small abs');
plot(t, M_v, 'LineWidth',1.4, 'DisplayName','medium abs');
plot(t, L_v, 'LineWidth',1.4, 'DisplayName','large abs');
plot(t, T_v, 'k', 'LineWidth',1.6, 'DisplayName','total abs');
xline(t_jump, 'k--', 'LineWidth',1.2);
grid on;
xlabel('Time (days)'); ylabel('Export (absolute)');
title('Legacy absolute export by class (volume)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'legacy_jump_fractions_and_absolute_export.png'));

fig2 = figure('Color','w');
subplot(2,1,1); hold on;
for b = bins05
    plot(t_win, Fv(win, b), 'LineWidth',1.2, 'DisplayName',sprintf('bin %d (D=%.3f cm)', b, d_cm_v(b)));
end
xline(t_jump, 'k--', 'LineWidth',1.2);
grid on;
xlabel('Time (days)'); ylabel('Flux (bin)');
title('Bins near 0.05 cm boundary');
legend('Location','best');

subplot(2,1,2); hold on;
for b = bins20
    plot(t_win, Fv(win, b), 'LineWidth',1.2, 'DisplayName',sprintf('bin %d (D=%.3f cm)', b, d_cm_v(b)));
end
xline(t_jump, 'k--', 'LineWidth',1.2);
grid on;
xlabel('Time (days)'); ylabel('Flux (bin)');
title('Bins near 0.20 cm boundary');
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'legacy_jump_boundary_bin_flux.png'));

mean_d_cm = sum(Fv .* d_cm_v', 2) ./ max(sum(Fv,2), eps);
fig3 = figure('Color','w');
plot(t, mean_d_cm*10, 'LineWidth',1.6); hold on;
xline(t_jump, 'k--', 'LineWidth',1.2);
grid on;
xlabel('Time (days)'); ylabel('Flux-weighted mean diameter (mm)');
title('Legacy flux-weighted mean diameter');
saveas(fig3, fullfile(fig_dir, 'legacy_jump_flux_weighted_mean_diameter.png'));

fprintf('\nSaved jump diagnostic figures to: %s\n', fig_dir);

