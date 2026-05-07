% run_may06_c3_sweep_matrix.m
% 4-case process matrix swept over c3 fragmentation strength.
%
% Cases: sink_only | sink+coag | sink+frag | sink+coag+frag
% Checks per case: neg_count, inventory_change_pct
%
% Goal: identify c3 where frag cases stay in a physically reasonable
% inventory change range (not +800%), while sink_only stays stable.
%
% Outputs:
%   Printed table + CSV in output/tables/
%   Figure in output/figures/

clear; close all; clc;
setup_paths

% ---- base config (match step-7 matrix settings) ----
base = SimulationConfig();
base.t_final    = 60;
base.delta_t    = 0.5;
base.n_sections = 20;
base.sinking_law   = 'kriest_8';
base.ds_kernel_mode = 'sinking_law';
base.box_depth  = 2000;   % m
base.r_to_rg    = 1.6;
base.enable_pp  = false;
base.pp_source  = 0;
base.growth     = 0;
base.enable_linear = false;
base.c4         = 1.45;
base.disagg_mode = 'legacy';

% ---- c3 sweep ----
c3_vals = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02];
n_c3    = numel(c3_vals);

% ---- case definitions ----
case_names = {'sink_only', 'sink_coag', 'sink_frag', 'sink_coag_frag'};
case_coag  = [false, true,  false, true];
case_frag  = [false, false, true,  true];
n_cases    = numel(case_names);

% ---- storage ----
neg_count_mat    = nan(n_c3, n_cases);
inv_change_mat   = nan(n_c3, n_cases);
converged_mat    = false(n_c3, n_cases);

% ---- header ----
fprintf('c3_sweep_matrix: t_final=%g days, n_sections=%d, kriest_8, box_depth=%g m\n', ...
    base.t_final, base.n_sections, base.box_depth);
fprintf('%-10s | %-16s | %-8s | %-18s\n', ...
    'c3', 'case', 'neg_cnt', 'inv_change_pct');
fprintf('%s\n', repmat('-', 1, 60));

% ---- run ----
for ic = 1:n_c3
    c3 = c3_vals(ic);
    for jc = 1:n_cases
        cfg = base.copy();
        cfg.c3             = c3;
        cfg.enable_coag    = case_coag(jc);
        cfg.enable_sinking = true;
        cfg.enable_disagg  = case_frag(jc);

        try
            sim = CoagulationSimulation(cfg);
            res = sim.run();

            Y = res.concentrations;

            neg_cnt = sum(Y(:) < -1e-30);

            N0 = sum(Y(1,:));
            Nf = sum(Y(end,:));
            if abs(N0) > 1e-60
                inv_pct = (Nf - N0) / N0 * 100;
            else
                inv_pct = NaN;
            end

            neg_count_mat(ic, jc)  = neg_cnt;
            inv_change_mat(ic, jc) = inv_pct;
            converged_mat(ic, jc)  = true;

            fprintf('%-10.4g | %-16s | %-8d | %-18.3f\n', ...
                c3, case_names{jc}, neg_cnt, inv_pct);

        catch ME
            fprintf('%-10.4g | %-16s | ERROR: %s\n', c3, case_names{jc}, ME.message);
        end
    end
end

% ---- recommendation ----
fprintf('\n=== Recommendation ===\n');
% For sink_frag case: find first c3 where inv_change < threshold
thresh = 50;  % percent; above this, fragmentation is dominating
frag_col = find(strcmp(case_names, 'sink_frag'));
if ~isempty(frag_col)
    ok_mask = inv_change_mat(:, frag_col) < thresh & converged_mat(:, frag_col);
    if any(ok_mask)
        best_idx = find(ok_mask, 1, 'last');
        fprintf('Largest c3 where sink_frag inv_change < %g%%: c3 = %.4g (inv_change = %.2f%%)\n', ...
            thresh, c3_vals(best_idx), inv_change_mat(best_idx, frag_col));
    else
        fprintf('No c3 in sweep gives sink_frag inv_change < %g%%. Extend sweep.\n', thresh);
    end
end

% ---- save CSV ----
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fileparts(script_dir);
tbl_dir    = fullfile(proj_root, 'output', 'tables');
if ~exist(tbl_dir, 'dir'), mkdir(tbl_dir); end

csv_path = fullfile(tbl_dir, 'may06_c3_sweep_matrix.csv');
fid = fopen(csv_path, 'w');
% header
hdr = 'c3';
for jc = 1:n_cases
    hdr = [hdr, sprintf(',neg_cnt_%s,inv_change_%s', case_names{jc}, case_names{jc})]; %#ok<AGROW>
end
fprintf(fid, '%s\n', hdr);
% rows
for ic = 1:n_c3
    row = sprintf('%.6g', c3_vals(ic));
    for jc = 1:n_cases
        row = [row, sprintf(',%.0f,%.4f', neg_count_mat(ic,jc), inv_change_mat(ic,jc))]; %#ok<AGROW>
    end
    fprintf(fid, '%s\n', row);
end
fclose(fid);
fprintf('\nCSV saved: %s\n', csv_path);

% ---- figure ----
fig_dir = fullfile(proj_root, 'output', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
fig = figure('Name', 'c3 sweep matrix: inventory change', 'Position', [100 100 800 500]);

for jc = 1:n_cases
    valid = converged_mat(:, jc) & isfinite(inv_change_mat(:, jc));
    if any(valid)
        semilogx(c3_vals(valid), inv_change_mat(valid, jc), 'o-', ...
            'LineWidth', 1.8, 'Color', colors(jc,:), ...
            'DisplayName', strrep(case_names{jc}, '_', '+'));
        hold on;
    end
end

yline(0,   'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(50,  'r:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(-50, 'b:',  'LineWidth', 1.2, 'HandleVisibility', 'off');

grid on;
xlabel('c3 (fragmentation strength)', 'FontSize', 12);
ylabel('Inventory change at t=60d (%)', 'FontSize', 12);
title(sprintf('4-case matrix: inventory change vs c3  (c4=%.2f, kriest8, H=%gm, t=%gd)', ...
    base.c4, base.box_depth, base.t_final), 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 10);

fname_fig = fullfile(fig_dir, 'may06_c3_sweep_matrix.png');
saveas(fig, fname_fig);
fprintf('Figure saved: %s\n', fname_fig);
