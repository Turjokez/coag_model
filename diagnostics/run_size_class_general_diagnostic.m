% run_size_class_general_diagnostic.m
% Build coarse size classes from the fine sectional model.
% This is an offline post-process check. No new ODE is solved here.

clear; close all; clc;
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
addpath(proj_root);
run(fullfile(proj_root, 'setup_paths.m'));
fig_dir = fullfile(proj_root, 'output', 'figures');
docs_fig_dir = fullfile(proj_root, 'docs', 'figures');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(docs_fig_dir, 'dir')
    mkdir(docs_fig_dir);
end

n_list = [2, 3, 4];
frag_eps = 1e-6;
clean_tol = 1e-10;

base = SimulationConfig();
base.t_final = 200;
base.delta_t = 1.0;
base.n_sections = 35;
base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;
base.enable_coag = true;
base.enable_sinking = true;
base.enable_linear = true;
base.growth = 0;
base.box_depth = 2000;
base.r_to_rg = 1.6;
base.sinking_law = 'kriest_8';
base.sinking_size = 'image';
base.sinking_scale = 0.8;

fprintf('\n=== No-frag test ===\n');
cfg_no = base.copy();
cfg_no.enable_disagg = false;
cfg_no.disagg_mode = 'legacy';

sim_no = CoagulationSimulation(cfg_no);
res_no = sim_no.run();

Y_no = res_no.concentrations;
v_lower = sim_no.grid.v_lower;
diam_i = res_no.output_data.diam_i(:);
betas = res_no.betas;
t_no = res_no.time(:);

coarse_no = struct([]);
for i = 1:numel(n_list)
    coarse_no(i).N = n_list(i);
    coarse_no(i).out = size_class_general(Y_no, v_lower, diam_i, betas, n_list(i));
    coarse_no(i).mass_err_max = max(abs(coarse_no(i).out.cons_mass_err));
    coarse_no(i).rate_err_max = max(abs(coarse_no(i).out.cons_rate_err));
    fprintf('N=%d | max |mass err| = %.3e | max |rate err| = %.3e\n', ...
        n_list(i), coarse_no(i).mass_err_max, coarse_no(i).rate_err_max);
end

exact_out = size_class_general(Y_no, v_lower, diam_i, betas, numel(v_lower), 1:(numel(v_lower)+1));
exact_state_err = max(abs(exact_out.C(:) - Y_no(:)));
exact_gain_err = max(abs(exact_out.gain(:) - exact_out.fine_gain(:)));
exact_loss_err = max(abs(exact_out.loss(:) - exact_out.fine_loss(:)));
exact_rate_err = max(abs(exact_out.cons_rate_err));

fprintf('Exact check N=n_sections | state err = %.3e | gain err = %.3e | loss err = %.3e | rate err = %.3e\n', ...
    exact_state_err, exact_gain_err, exact_loss_err, exact_rate_err);

frag_ok = all([coarse_no.mass_err_max] < clean_tol) && all([coarse_no.rate_err_max] < clean_tol);
coarse_frag = struct([]);
res_frag = [];
t_frag = [];

if frag_ok
    fprintf('\n=== Frag test (eps=%0.0e) ===\n', frag_eps);
    cfg_frag = base.copy();
    cfg_frag.enable_disagg = true;
    cfg_frag.disagg_mode = 'operator_split';
    cfg_frag.disagg_epsilon = frag_eps;
    cfg_frag.disagg_outer_dt = 1/24;

    sim_frag = CoagulationSimulation(cfg_frag);
    res_frag = sim_frag.run();
    Y_frag = res_frag.concentrations;
    t_frag = res_frag.time(:);

    for i = 1:numel(n_list)
        coarse_frag(i).N = n_list(i);
        coarse_frag(i).out = size_class_general(Y_frag, v_lower, diam_i, betas, n_list(i));
        coarse_frag(i).mass_err_max = max(abs(coarse_frag(i).out.cons_mass_err));
        coarse_frag(i).rate_err_max = max(abs(coarse_frag(i).out.cons_rate_err));
        fprintf('Frag N=%d | max |mass err| = %.3e | max |rate err| = %.3e\n', ...
            n_list(i), coarse_frag(i).mass_err_max, coarse_frag(i).rate_err_max);
    end
else
    fprintf('\nNo-frag conservation was not clean enough, so frag check was skipped.\n');
end

failure_log = run_failure_checks(Y_no, v_lower, diam_i, betas);

fig1 = figure('Color', 'w', 'Position', [80 80 1400 420]);
D = res_no.output_data.diam_i(:);
Nf = res_no.output_data.nspec_i(end, :)';
for i = 1:numel(n_list)
    subplot(1, numel(n_list), i);
    plot_grouped_spectrum(D, Nf, coarse_no(i).out, sprintf('Final fine spectrum with N=%d groups', n_list(i)));
end
saveas(fig1, fullfile(fig_dir, 'size_class_general_spectrum_groups.png'));

fig2 = figure('Color', 'w', 'Position', [80 80 1400 420]);
for i = 1:numel(n_list)
    subplot(1, numel(n_list), i);
    plot_coarse_amounts(t_no, coarse_no(i).out, sprintf('Coarse amount over time (N=%d)', n_list(i)));
end
saveas(fig2, fullfile(fig_dir, 'size_class_general_amount_timeseries.png'));

fig3 = figure('Color', 'w', 'Position', [80 80 1500 760]);
pick_no = 1;
subplot(2, 3, 1);
plot_terms(t_no, coarse_no(pick_no).out.gain, sprintf('No-frag gain (N=%d)', n_list(pick_no)));
subplot(2, 3, 2);
plot_terms(t_no, coarse_no(pick_no).out.loss, sprintf('No-frag loss (N=%d)', n_list(pick_no)));
subplot(2, 3, 3);
plot_terms(t_no, coarse_no(pick_no).out.net, sprintf('No-frag net dC/dt (N=%d)', n_list(pick_no)));
if frag_ok
    pick_frag = 1;
    subplot(2, 3, 4);
    plot_terms(t_frag, coarse_frag(pick_frag).out.gain, sprintf('Frag gain (N=%d)', n_list(pick_frag)));
    subplot(2, 3, 5);
    plot_terms(t_frag, coarse_frag(pick_frag).out.loss, sprintf('Frag loss (N=%d)', n_list(pick_frag)));
    subplot(2, 3, 6);
    plot_terms(t_frag, coarse_frag(pick_frag).out.net, sprintf('Frag net dC/dt (N=%d)', n_list(pick_frag)));
else
    pick_hi = numel(n_list);
    subplot(2, 3, 4);
    plot_terms(t_no, coarse_no(pick_hi).out.gain, sprintf('No-frag gain (N=%d)', n_list(pick_hi)));
    subplot(2, 3, 5);
    plot_terms(t_no, coarse_no(pick_hi).out.loss, sprintf('No-frag loss (N=%d)', n_list(pick_hi)));
    subplot(2, 3, 6);
    plot_terms(t_no, coarse_no(pick_hi).out.net, sprintf('No-frag net dC/dt (N=%d)', n_list(pick_hi)));
end
saveas(fig3, fullfile(fig_dir, 'size_class_general_gain_loss_net.png'));

fig4 = figure('Color', 'w', 'Position', [80 80 1350 420]);
subplot(1, 3, 1);
plot_error_summary(n_list, coarse_no, coarse_frag, 'mass_err_max', frag_ok, 'Max amount error');
subplot(1, 3, 2);
plot_error_summary(n_list, coarse_no, coarse_frag, 'rate_err_max', frag_ok, 'Max coag-net error');
subplot(1, 3, 3);
plot_summary_text(exact_state_err, exact_gain_err, exact_loss_err, exact_rate_err, failure_log, frag_ok, frag_eps);
saveas(fig4, fullfile(fig_dir, 'size_class_general_conservation.png'));

copy_one(fullfile(fig_dir, 'size_class_general_spectrum_groups.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'size_class_general_amount_timeseries.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'size_class_general_gain_loss_net.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'size_class_general_conservation.png'), docs_fig_dir);

summary = struct();
summary.no_frag = coarse_no;
summary.frag_ok = frag_ok;
summary.frag = coarse_frag;
summary.exact_state_err = exact_state_err;
summary.exact_gain_err = exact_gain_err;
summary.exact_loss_err = exact_loss_err;
summary.exact_rate_err = exact_rate_err;
summary.failure_log = failure_log;
summary.fig_dir = fig_dir;
summary.docs_fig_dir = docs_fig_dir;

save(fullfile(fig_dir, 'size_class_general_summary.mat'), 'summary');

fprintf('\nSaved figures to: %s\n', fig_dir);
fprintf('Copied report figures to: %s\n', docs_fig_dir);

function plot_grouped_spectrum(D, Nf, out, ttl)
cols = lines(size(out.groups, 1));
loglog(D, Nf, 'k-', 'LineWidth', 1.5); hold on;
yl = [max(min(Nf(Nf > 0)), eps), max(Nf)];
for ic = 1:size(out.groups, 1)
    idx = out.groups(ic, :);
    d_lo = D(idx(1));
    d_hi = D(idx(2));
    patch([d_lo d_hi d_hi d_lo], [yl(1) yl(1) yl(2) yl(2)], cols(ic,:), ...
        'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(out.diam_mid(ic), yl(2) / (1.4 + 0.15 * ic), sprintf('C%d', ic), ...
        'HorizontalAlignment', 'center', 'Color', cols(ic,:));
end
for ic = 1:(size(out.groups, 1) - 1)
    i_right = out.groups(ic, 2);
    x_sep = sqrt(D(i_right) * D(i_right + 1));
    xline(x_sep, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Image diameter (cm)');
ylabel('Number spectrum N');
title(ttl);
xlim([min(D) max(D)]);
end

function plot_coarse_amounts(t, out, ttl)
cols = lines(size(out.C, 2));
for ic = 1:size(out.C, 2)
    plot(t, out.C(:, ic), 'LineWidth', 1.3, 'Color', cols(ic, :), ...
        'DisplayName', sprintf('class %d', ic));
    hold on;
end
grid on;
xlabel('Time (days)');
ylabel('Coarse class amount C');
title(ttl);
legend('Location', 'best');
end

function plot_terms(t, A, ttl)
cols = lines(size(A, 2));
for ic = 1:size(A, 2)
    plot(t, A(:, ic), 'LineWidth', 1.2, 'Color', cols(ic, :), ...
        'DisplayName', sprintf('class %d', ic));
    hold on;
end
grid on;
xlabel('Time (days)');
ylabel('Rate');
title(ttl);
legend('Location', 'best');
end

function plot_error_summary(n_list, coarse_no, coarse_frag, field_name, frag_ok, ttl)
x = 1:numel(n_list);
y_no = zeros(size(x));
for i = 1:numel(n_list)
    y_no(i) = coarse_no(i).(field_name);
end
semilogy(x, y_no, '-o', 'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', 'no-frag');
hold on;
if frag_ok
    y_fg = zeros(size(x));
    for i = 1:numel(n_list)
        y_fg(i) = coarse_frag(i).(field_name);
    end
    semilogy(x, y_fg, '-s', 'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', 'frag');
end
grid on;
xlabel('Coarse class count N');
ylabel('Max absolute error');
title(ttl);
legend('Location', 'best');
xticks(x);
xticklabels(string(n_list));
ylim([1e-22 1e-18]);
end

function plot_summary_text(exact_state_err, exact_gain_err, exact_loss_err, exact_rate_err, failure_log, frag_ok, frag_eps)
axis off;
txt = {
    'Exact check'
    sprintf('state err = %.2e', exact_state_err)
    sprintf('gain err = %.2e', exact_gain_err)
    sprintf('loss err = %.2e', exact_loss_err)
    sprintf('rate err = %.2e', exact_rate_err)
    ' '
    'Failure checks'
    sprintf('bad N < 2 = %d', failure_log.bad_N)
    sprintf('repeat edges = %d', failure_log.bad_edges_repeat)
    sprintf('short edges = %d', failure_log.bad_edges_short)
    sprintf('non-monotone v_{lower} = %d', failure_log.bad_v_lower)
    ' '
    sprintf('frag test used eps = %.0e', frag_eps)
    sprintf('frag run included = %d', frag_ok)
    };
text(0.02, 0.98, txt, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontName', 'Helvetica', 'FontSize', 11);
title('Check summary');
end

function log_out = run_failure_checks(Y, v_lower, diam_i, betas)
n_bad = numel(v_lower);
log_out = struct();
log_out.bad_N = did_fail(@() size_class_general(Y, v_lower, diam_i, betas, 1));
log_out.bad_edges_repeat = did_fail(@() size_class_general(Y, v_lower, diam_i, betas, 3, [1 10 10 (n_bad + 1)]));
log_out.bad_edges_short = did_fail(@() size_class_general(Y, v_lower, diam_i, betas, 3, [1 10 (n_bad + 1)]));
bad_v = v_lower(1:6);
bad_v(end) = bad_v(end-1);
log_out.bad_v_lower = did_fail(@() size_class_general(Y(:, 1:6), bad_v, diam_i(1:6), betas, 2));

fprintf('\nFailure checks:\n');
fprintf('bad N < 2: %d\n', log_out.bad_N);
fprintf('repeat edges: %d\n', log_out.bad_edges_repeat);
fprintf('short edges: %d\n', log_out.bad_edges_short);
fprintf('non-monotone v_lower: %d\n', log_out.bad_v_lower);
end

function tf = did_fail(fh)
tf = false;
try
    fh();
catch
    tf = true;
end
end

function copy_one(src_file, dst_dir)
[~, name, ext] = fileparts(src_file);
copyfile(src_file, fullfile(dst_dir, [name, ext]));
end
