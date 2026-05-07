% run_may06_operator_split_matrix.m
% 4-case process matrix using operator_split disaggregation.
% Sweeps D_max [cm] and compares against legacy c3 result.
%
% Operator_split is mass-conserving: particles above D_max are
% redistributed to smaller bins, total biovolume is preserved.
% This is the physically preferred mode for the 1-D depth model.
%
% Checks: neg_count, inventory_change_pct, biovolume_change_pct
% (biovolume should be conserved by operator_split disagg).
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
base.sinking_law    = 'kriest_8';
base.ds_kernel_mode = 'sinking_law';
base.box_depth  = 2000;   % m
base.r_to_rg    = 1.6;
base.enable_pp  = false;
base.pp_source  = 0;
base.growth     = 0;
base.enable_linear = false;

% operator_split settings
base.disagg_mode      = 'operator_split';
base.disagg_outer_dt  = 1/24;   % 1 hour
base.disagg_frac_next = 2/3;
base.disagg_epsilon   = [];      % use dmax_cm override

% ---- D_max sweep [cm] ----
dmax_vals = [0.5, 1.0, 2.0, 5.0];   % cm  (0.5cm=5mm to 5cm=50mm)
n_dmax    = numel(dmax_vals);

% ---- case definitions ----
case_names = {'sink_only', 'sink_coag', 'sink_frag', 'sink_coag_frag'};
case_coag  = [false, true,  false, true];
case_frag  = [false, false, true,  true];
n_cases    = numel(case_names);

% ---- storage ----
neg_count_arr    = nan(n_dmax, n_cases);
inv_change_arr   = nan(n_dmax, n_cases);
biovol_change_arr = nan(n_dmax, n_cases);
converged_arr    = false(n_dmax, n_cases);

% ---- header ----
fprintf('operator_split_matrix: t_final=%g d, n_sec=%d, kriest_8, H=%g m\n', ...
    base.t_final, base.n_sections, base.box_depth);
fprintf('%-8s | %-16s | %-8s | %-20s | %-22s\n', ...
    'dmax_cm', 'case', 'neg_cnt', 'inv_change_pct', 'biovol_change_pct');
fprintf('%s\n', repmat('-', 1, 80));

% ---- run ----
for id = 1:n_dmax
    dmax_cm = dmax_vals(id);
    for jc = 1:n_cases
        cfg = base.copy();
        cfg.disagg_dmax_cm = dmax_cm;
        cfg.enable_coag    = case_coag(jc);
        cfg.enable_sinking = true;
        cfg.enable_disagg  = case_frag(jc);

        try
            sim = CoagulationSimulation(cfg);
            res = sim.run();

            Y      = res.concentrations;
            av_vol = sim.grid.av_vol(:)';

            neg_cnt = sum(Y(:) < -1e-30);

            % inventory (sum of state = number-like quantity)
            N0 = sum(Y(1,:));
            Nf = sum(Y(end,:));
            if abs(N0) > 1e-60
                inv_pct = (Nf - N0) / N0 * 100;
            else
                inv_pct = NaN;
            end

            % biovolume (conserved by operator_split)
            B0 = Y(1,:)   * av_vol';
            Bf = Y(end,:) * av_vol';
            if abs(B0) > 1e-60
                biovol_pct = (Bf - B0) / B0 * 100;
            else
                biovol_pct = NaN;
            end

            neg_count_arr(id, jc)     = neg_cnt;
            inv_change_arr(id, jc)    = inv_pct;
            biovol_change_arr(id, jc) = biovol_pct;
            converged_arr(id, jc)     = true;

            fprintf('%-8.2f | %-16s | %-8d | %-20.3f | %-22.6f\n', ...
                dmax_cm, case_names{jc}, neg_cnt, inv_pct, biovol_pct);

        catch ME
            fprintf('%-8.2f | %-16s | ERROR: %s\n', dmax_cm, case_names{jc}, ME.message);
        end
    end
end

% ---- summary ----
fprintf('\n=== operator_split verification ===\n');
fprintf('For frag cases, biovolume_change should be ~same as sink_only\n');
fprintf('(operator_split disagg conserves biovolume, only sinking removes it)\n\n');

% Check: for sink_frag, biovol_change should track sink_only closely
sink_col = find(strcmp(case_names, 'sink_only'));
frag_col = find(strcmp(case_names, 'sink_frag'));
if ~isempty(sink_col) && ~isempty(frag_col)
    fprintf('%-8s | %-22s | %-22s | %-16s\n', ...
        'dmax_cm', 'sink_only biovol%', 'sink_frag biovol%', 'difference');
    fprintf('%s\n', repmat('-', 1, 72));
    for id = 1:n_dmax
        bs = biovol_change_arr(id, sink_col);
        bf = biovol_change_arr(id, frag_col);
        diff_pct = abs(bf - bs);
        fprintf('%-8.2f | %-22.4f | %-22.4f | %-16.6f\n', ...
            dmax_vals(id), bs, bf, diff_pct);
    end
end

% ---- save CSV ----
script_dir = fileparts(mfilename('fullpath'));
proj_root  = fileparts(script_dir);
tbl_dir    = fullfile(proj_root, 'output', 'tables');
if ~exist(tbl_dir, 'dir'), mkdir(tbl_dir); end

csv_path = fullfile(tbl_dir, 'may06_operator_split_matrix.csv');
fid = fopen(csv_path, 'w');
hdr = 'dmax_cm';
for jc = 1:n_cases
    hdr = [hdr, sprintf(',neg_cnt_%s,inv_change_%s,biovol_change_%s', ...
        case_names{jc}, case_names{jc}, case_names{jc})]; %#ok<AGROW>
end
fprintf(fid, '%s\n', hdr);
for id = 1:n_dmax
    row = sprintf('%.4g', dmax_vals(id));
    for jc = 1:n_cases
        row = [row, sprintf(',%.0f,%.4f,%.6f', ...
            neg_count_arr(id,jc), inv_change_arr(id,jc), biovol_change_arr(id,jc))]; %#ok<AGROW>
    end
    fprintf(fid, '%s\n', row);
end
fclose(fid);
fprintf('\nCSV saved: %s\n', csv_path);

% ---- figure: biovol_change for frag cases by dmax ----
fig_dir = fullfile(proj_root, 'output', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
fig = figure('Name', 'Operator-split matrix: biovolume change', 'Position', [100 100 800 500]);

for jc = 1:n_cases
    valid = converged_arr(:, jc) & isfinite(biovol_change_arr(:, jc));
    if any(valid)
        plot(dmax_vals(valid), biovol_change_arr(valid, jc), 'o-', ...
            'LineWidth', 1.8, 'Color', colors(jc,:), ...
            'DisplayName', strrep(case_names{jc}, '_', '+'));
        hold on;
    end
end

yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on;
xlabel('D_{max} [cm]', 'FontSize', 12);
ylabel('Biovolume change at t=60d (%)', 'FontSize', 12);
title(sprintf('Operator-split matrix: biovolume change vs D_{max}  (kriest8, H=%gm, t=%gd)', ...
    base.box_depth, base.t_final), 'FontSize', 11);
legend('Location', 'best', 'FontSize', 10);

fname_fig = fullfile(fig_dir, 'may06_operator_split_matrix.png');
saveas(fig, fname_fig);
fprintf('Figure saved: %s\n', fname_fig);
