% run_may06_frag_rate_analysis.m
% Analytical comparison of fragmentation rate vs sinking rate per bin.
% No simulation needed. Uses grid math only.
%
% Purpose: show WHY c3=0.02 causes +869% number change.
% For the legacy disagg term:  rate(isec) = c3 * c4^isec  [1/day]
% For sinking:                 rate(isec) = w_i / H        [1/day]
% When frag_rate >> sink_rate, fragmentation dominates and number explodes.
%
% Outputs:
%   Printed table per bin (frag rate, sink rate, ratio) for each c3
%   Figure saved to output/figures/

clear; close all; clc;
setup_paths

% ---- config (match step-7 matrix settings) ----
cfg = SimulationConfig();
cfg.n_sections  = 20;
cfg.sinking_law = 'kriest_8';
cfg.box_depth   = 2000;   % m
cfg.r_to_rg     = 1.6;
cfg.c4          = 1.45;

grid = DerivedGrid(cfg);

% ---- sinking rate per bin [1/day] ----
v_cms  = SettlingVelocityService.velocityForSections(grid, cfg);  % cm/s
w_mday = (v_cms / 100) * cfg.day_to_sec;                          % m/day
sink_rate = w_mday / cfg.box_depth;                               % 1/day

% ---- c3 values to compare ----
c3_vals = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02];
c4      = cfg.c4;
isec    = (1:cfg.n_sections)';

% ---- diameters for reference ----
d_mm = grid.getVolumeDiameters() * 10;  % cm -> mm

% ---- print table for current default c3=0.02 ----
fprintf('\n=== Rate comparison: legacy disagg vs sinking (c3=0.02, c4=%.2f, H=%g m) ===\n', ...
    c4, cfg.box_depth);
fprintf('%-4s  %-8s  %-12s  %-12s  %-10s\n', ...
    'bin', 'd_mm', 'frag_rate/d', 'sink_rate/d', 'ratio');
fprintf('%s\n', repmat('-', 1, 55));

c3_default = 0.02;
for i = 1:cfg.n_sections
    fr = c3_default * c4^i;
    sr = sink_rate(i);
    ratio = fr / max(sr, 1e-30);
    flag = '';
    if ratio > 10
        flag = ' <<< FRAGMENTATION DOMINATES';
    elseif ratio > 1
        flag = ' < frag > sink';
    end
    fprintf('%-4d  %-8.3f  %-12.4e  %-12.4e  %-10.2f%s\n', ...
        i, d_mm(i), fr, sr, ratio, flag);
end

% ---- find c3 where max ratio across all bins is ~1 (balanced) ----
fprintf('\n=== c3 balance check: max(frag_rate / sink_rate) across all bins ===\n');
fprintf('%-10s  %-16s  %-8s\n', 'c3', 'max_ratio', 'balanced_bin');
fprintf('%s\n', repmat('-', 1, 40));

max_ratios = zeros(size(c3_vals));
for ic = 1:numel(c3_vals)
    c3  = c3_vals(ic);
    fr  = c3 * c4.^isec;
    rat = fr ./ max(sink_rate(:), 1e-30);
    max_ratios(ic) = max(rat);
    [~, imax] = max(rat);
    fprintf('%-10.4g  %-16.3e  bin %-2d (d=%.3f mm)\n', ...
        c3, max_ratios(ic), imax, d_mm(imax));
end

% NOTE: max_ratio is always at the largest bin (bin 20) where c4^20 is large.
% But bin 20 has near-zero population in the initial spectrum (10^{-19} relative),
% so the actual fragmentation flux there is negligible.
% The physically meaningful check is the ratio at a representative mid-range bin
% where population is significant (bins 5-10 for a typical decaying spectrum).

% ---- ratio at representative bin 8 (d~0.11mm, significant population) ----
rep_bin = 8;
fprintf('\n=== Representative-bin check (bin %d, d=%.3f mm) ===\n', rep_bin, d_mm(rep_bin));
fprintf('%-10s  %-16s  %-10s\n', 'c3', 'frag_rate/d', 'frag/sink');
fprintf('%s\n', repmat('-', 1, 40));

rep_sink = sink_rate(rep_bin);
rep_ratios = zeros(size(c3_vals));
for ic = 1:numel(c3_vals)
    c3 = c3_vals(ic);
    fr = c3 * c4^rep_bin;
    rep_ratios(ic) = fr / max(rep_sink, 1e-30);
    fprintf('%-10.4g  %-16.4e  %-10.2f\n', c3, fr, rep_ratios(ic));
end

% Recommended c3: ratio at representative bin < 2 (frag weaker than 2x sinking)
safe_rep = rep_ratios < 2;
if any(safe_rep)
    rec_idx = find(safe_rep, 1, 'last');
    fprintf('\nRecommended c3 (frag/sink < 2 at bin %d): c3 = %.4g\n', rep_bin, c3_vals(rec_idx));
    fprintf('  -> Run run_may06_c3_sweep_matrix.m to confirm with actual simulation.\n');
else
    fprintf('\nNo c3 in sweep gives frag/sink < 2 at bin %d. Extend sweep.\n', rep_bin);
end

% ---- figure: frag rate and sink rate per bin for several c3 values ----
plot_c3 = [0.0002, 0.001, 0.005, 0.02];
colors   = lines(numel(plot_c3));

fig = figure('Name', 'Frag vs Sink rate per bin', 'Position', [100 100 900 500]);

semilogy(isec, sink_rate, 'k-', 'LineWidth', 2.5, 'DisplayName', 'sinking rate (w/H)');
hold on;
for ic = 1:numel(plot_c3)
    c3 = plot_c3(ic);
    fr = c3 * c4.^isec;
    semilogy(isec, fr, '--', 'LineWidth', 1.8, 'Color', colors(ic,:), ...
        'DisplayName', sprintf('frag rate c3=%.4g', c3));
end

grid on;
xlabel('Section (bin index)', 'FontSize', 12);
ylabel('Rate coefficient [1/day]', 'FontSize', 12);
title(sprintf('Fragmentation vs Sinking rate per bin  (kriest8, H=%g m)', cfg.box_depth), ...
    'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 10);
xlim([1, cfg.n_sections]);

% Add secondary x-axis labels for diameter
xt = [1, 5, 10, 15, 20];
xt = xt(xt <= cfg.n_sections);
ax = gca;
ax.XTick = xt;
ax.XTickLabel = arrayfun(@(i) sprintf('%d\n(%.2fmm)', i, d_mm(i)), xt, 'UniformOutput', false);

% ---- save figure ----
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fileparts(script_dir);
outdir     = fullfile(proj_root, 'output', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

fname = fullfile(outdir, 'may06_frag_rate_analysis.png');
saveas(fig, fname);
fprintf('\nFigure saved: %s\n', fname);
