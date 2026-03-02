% run_day63_numeric_check.m
% Check if day~63 jump is numeric issue or model behavior.

clear; close all; clc;
setup_paths

% ---------- baseline legacy config ----------
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
cfg.c4 = 1.45;
cfg.r_to_rg = 1.6;
cfg.sinking_law = 'kriest_8';
cfg.box_depth = 2000;
cfg.disagg_mode = 'legacy';

fprintf('Running baseline legacy case...\n');
sim0 = CoagulationSimulation(cfg);
res0 = sim0.run();
out0 = analyze_jump(res0, cfg, 'baseline');

% A looser solver test: if jump is pure solver noise, this should move a lot.
fprintf('Running loose solver test...\n');
sim1 = CoagulationSimulation(cfg);
loose_opt = struct('RelTol', 1e-8, 'AbsTol', 1e-12);
res1 = sim1.run('solver_options', loose_opt);
out1 = analyze_jump(res1, cfg, 'loose');

fprintf('\n=== Day-63 Numeric Check Summary ===\n');
print_summary(out0);
print_summary(out1);

dt_jump = abs(out0.t_jump - out1.t_jump);
dmetric = abs(out0.jump_metric - out1.jump_metric);
fprintf('\nJump-time change (baseline vs loose): %.4f days\n', dt_jump);
fprintf('Jump-metric change (baseline vs loose): %.6f\n', dmetric);

% ---------- save one quick plot ----------
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
fig_dir = fullfile(proj_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

fig = figure('Color','w');
subplot(2,1,1); hold on;
plot(out0.t, out0.fS, 'LineWidth',1.3, 'DisplayName','small (baseline)');
plot(out0.t, out0.fM, 'LineWidth',1.3, 'DisplayName','medium (baseline)');
plot(out0.t, out0.fL, 'LineWidth',1.3, 'DisplayName','large (baseline)');
xline(out0.t_jump, 'k--', 'DisplayName','jump');
grid on; xlabel('Time (days)'); ylabel('Export fraction');
title('Legacy fractions (baseline)');
legend('Location','best');

subplot(2,1,2); hold on;
plot(out0.t_window, out0.top_bins_window, 'LineWidth',1.2);
xline(out0.t_jump, 'k--');
set(gca, 'YScale', 'log');
grid on; xlabel('Time (days)'); ylabel('Concentration (log)');
title('Top 8 bin concentrations near jump (baseline)');

saveas(fig, fullfile(fig_dir, 'day63_numeric_check.png'));
fprintf('\nSaved figure: %s\n', fullfile(fig_dir, 'day63_numeric_check.png'));

% ---------- local helpers ----------
function out = analyze_jump(res, cfg, tag)
    t = res.time(:);
    Y = res.concentrations;
    Fv = max(res.output_data.fluxspec, 0);
    d_cm = res.output_data.diam_v(:);

    n_t = numel(t);
    fS = zeros(n_t,1); fM = fS; fL = fS; T = fS;
    for k = 1:n_t
        [fS(k), fM(k), fL(k), tot] = export_size_classes(d_cm, Fv(k,:).', 'volume');
        T(k) = tot.total;
    end

    df = abs(diff(fS)) + abs(diff(fM)) + abs(diff(fL));
    [jump_metric, ix] = max(df);
    t_jump = t(ix+1);

    y0 = Y(ix,:).';
    y1 = Y(ix+1,:).';
    n = numel(y0);

    % disagg sign term: v(i) - c4*v(i+1), i=2..n-1
    s0 = y0(2:n-1) - cfg.c4 * y0(3:n);
    s1 = y1(2:n-1) - cfg.c4 * y1(3:n);
    sign_flip_idx = find(sign(s0) ~= sign(s1)) + 1;

    top_idx = (n-7):n;
    t_win_mask = (t >= t_jump-2) & (t <= t_jump+2);

    out = struct();
    out.tag = tag;
    out.t = t;
    out.fS = fS; out.fM = fM; out.fL = fL;
    out.jump_metric = jump_metric;
    out.t_jump = t_jump;
    out.t0 = t(ix); out.t1 = t(ix+1);
    out.frac0 = [fS(ix), fM(ix), fL(ix)];
    out.frac1 = [fS(ix+1), fM(ix+1), fL(ix+1)];
    out.total0 = T(ix); out.total1 = T(ix+1);
    out.any_nan = any(~isfinite(Y(:))) || any(~isfinite(Fv(:)));
    out.min_pos = min(Y(Y > 0));
    out.max_val = max(Y(:));
    out.n_below_eps_pre = sum(abs(y0) < eps);
    out.n_below_eps_post = sum(abs(y1) < eps);
    out.sign_flip_idx = sign_flip_idx(:).';
    out.top_idx = top_idx;
    out.top_pre = y0(top_idx).';
    out.top_post = y1(top_idx).';
    out.t_window = t(t_win_mask);
    out.top_bins_window = Y(t_win_mask, top_idx);
end

function print_summary(o)
    fprintf('\n[%s]\n', o.tag);
    fprintf('jump metric = %.6f at %.2f -> %.2f days\n', o.jump_metric, o.t0, o.t1);
    fprintf('fractions before = [%.3f %.3f %.3f]\n', o.frac0(1), o.frac0(2), o.frac0(3));
    fprintf('fractions after  = [%.3f %.3f %.3f]\n', o.frac1(1), o.frac1(2), o.frac1(3));
    fprintf('total export before/after = %.3e -> %.3e\n', o.total0, o.total1);
    fprintf('any NaN/Inf = %d | min positive = %.3e | max val = %.3e\n', ...
        o.any_nan, o.min_pos, o.max_val);
    fprintf('count(|y|<eps) pre/post = %d / %d\n', o.n_below_eps_pre, o.n_below_eps_post);
    fprintf('sign flips in disagg term bins = ');
    if isempty(o.sign_flip_idx)
        fprintf('none\n');
    else
        fprintf('%s\n', mat2str(o.sign_flip_idx));
    end
    fprintf('top bins pre  = %s\n', mat2str(o.top_pre, 3));
    fprintf('top bins post = %s\n', mat2str(o.top_post, 3));
end
