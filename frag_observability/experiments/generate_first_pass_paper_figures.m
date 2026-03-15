function out = generate_first_pass_paper_figures(opts)
% generate_first_pass_paper_figures
% Build first-pass paper figures from the frag observability framework.

if nargin < 1 || isempty(opts)
    opts = struct();
end

exp_dir = fileparts(mfilename('fullpath'));
frag_root = fileparts(exp_dir);
repo_root = fileparts(frag_root);

addpath(exp_dir);
addpath(fullfile(frag_root, 'diagnostics'));
rehash;

opts = set_default(opts, 'apply_uvp', true);
opts = set_default(opts, 'uvp_sample_L', 1.0);
opts = set_default(opts, 'fit_range_cm', [0.01, 0.6]);
opts = set_default(opts, 'plot_range_cm', [0.005, 2.0]);
opts = set_default(opts, 'output_root', fullfile(repo_root, 'frag_observability'));
opts = set_default(opts, 'use_saved_results', true);
opts = set_default(opts, 'sample_compare_L', [1, 10]);
opts = set_default(opts, 'sample_repeats', 25);
opts = set_default(opts, 'random_seed', 1);

fig_dir = fullfile(opts.output_root, 'figures');
note_dir = fullfile(opts.output_root, 'notes');
tbl_dir = fullfile(opts.output_root, 'tables');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(note_dir, 'dir'), mkdir(note_dir); end

save_file = fullfile(tbl_dir, 'frag_observability_case_summary.mat');
if opts.use_saved_results && exist(save_file, 'file')
    tmp = load(save_file, 'all_runs', 'opts');
    all_runs = tmp.all_runs;
    res = struct();
    res.all_runs = all_runs;
    if isfield(tmp, 'opts')
        opts.fit_range_cm = tmp.opts.fit_range_cm;
        opts.plot_range_cm = tmp.opts.plot_range_cm;
        opts.apply_uvp = tmp.opts.apply_uvp;
        opts.uvp_sample_L = tmp.opts.uvp_sample_L;
    end
else
    res = run_frag_observability_framework(opts);
    all_runs = res.all_runs;
end

fit_range_cm = opts.fit_range_cm;
plot_range_cm = opts.plot_range_cm;

% ---------- Figure 1: (b, kappa) separability ----------
fig1 = figure('Color', 'w');
tl1 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot_b_kappa_separability(gather_fits(all_runs, 'fit_img'), gca, false);
title('Image diameter | color = sinking law | shape = frag level');

nexttile;
plot_b_kappa_separability(gather_fits(all_runs, 'fit_vol'), gca, false);
title('Volume diameter | color = sinking law | shape = frag level');

title(tl1, '(b, kappa) separability');
saveas(fig1, fullfile(fig_dir, 'fig01_b_kappa_separability.png'));

% ---------- Figure 2: baseline PSD + fit + residuals ----------
ref_law_idx = pick_reference_law(all_runs);
ref_fit = all_runs(ref_law_idx).cases(1).fit_img;

fig2 = figure('Color', 'w');
tl2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
loglog(ref_fit.Dall, ref_fit.Nall, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.7);
hold on;
loglog(ref_fit.Dsel, ref_fit.Nsel, 'o', 'Color', [0.15 0.15 0.15], ...
    'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 4, ...
    'DisplayName', 'fit bins');
loglog(ref_fit.Dfit, ref_fit.Nfit, '--', 'Color', [0.12 0.45 0.78], ...
    'LineWidth', 1.4, 'DisplayName', sprintf('power-law fit | b=%.3f', ref_fit.b));
grid on;
xlabel('Image diameter (cm)');
ylabel('N(D)');
title(sprintf('Baseline PSD | %s | no-frag', all_runs(ref_law_idx).law), 'Interpreter', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(plot_range_cm);
legend('Location', 'best');

nexttile;
plot(ref_fit.x, ref_fit.resid, '-o', 'Color', [0.12 0.45 0.78], ...
    'LineWidth', 1.4, 'MarkerSize', 4, ...
    'DisplayName', sprintf('residuals | kappa=%.3g', ref_fit.kappa));
hold on;
plot(ref_fit.x, ref_fit.resid_quad, '--', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.3, 'DisplayName', 'quadratic fit');
yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on;
xlabel('log10(D)');
ylabel('Residual');
title('Residual shape in fit window');
legend('Location', 'best');

title(tl2, 'Baseline PSD, fit, and residuals');
saveas(fig2, fullfile(fig_dir, 'fig02_baseline_psd_fit_residuals.png'));

% ---------- Figure 3: kappa vs fragmentation strength ----------
fig3 = figure('Color', 'w');
tl3 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot_kappa_vs_level(all_runs, 'fit_img');
xlabel('Fragmentation level');
ylabel('kappa');
title('Image diameter');
grid on;
legend('Location', 'best');

nexttile; hold on;
plot_kappa_vs_level(all_runs, 'fit_vol');
xlabel('Fragmentation level');
ylabel('kappa');
title('Volume diameter');
grid on;
legend('Location', 'best');

title(tl3, 'kappa vs fragmentation strength');
saveas(fig3, fullfile(fig_dir, 'fig03_kappa_vs_fragmentation_strength.png'));

% ---------- Figure 4: clean vs observed window vs noisy ----------
[law_idx, frag_idx] = strongest_case_idx(all_runs);
law_name = all_runs(law_idx).law;
case_clean0 = all_runs(law_idx).cases(1);
case_clean1 = all_runs(law_idx).cases(frag_idx);

fig4 = figure('Color', 'w');
tl4 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% panel 1: clean full PSD
nexttile; hold on;
plot_clean_pair(case_clean0.fit_img, case_clean1.fit_img, true);
title(sprintf('Clean PSD | %s', law_name), 'Interpreter', 'none');
xlim(plot_range_cm);

% panel 2: observed window only
nexttile; hold on;
plot_window_pair(case_clean0.fit_img, case_clean1.fit_img, fit_range_cm);
title('Observed window');

% panel 3: noisy sampled PSD
nexttile; hold on;
if opts.apply_uvp && isfield(case_clean0, 'fit_img_sampled') && case_clean0.fit_img_sampled.ok ...
        && isfield(case_clean1, 'fit_img_sampled') && case_clean1.fit_img_sampled.ok
    plot_noisy_pair(case_clean0.fit_img, case_clean0.fit_img_sampled, ...
        case_clean1.fit_img, case_clean1.fit_img_sampled, opts.uvp_sample_L);
else
    plot_clean_pair(case_clean0.fit_img, case_clean1.fit_img, false);
    title('Noisy view unavailable');
end

title(tl4, sprintf('Clean vs observed-window vs noisy | strongest case: %s %s', ...
    law_name, case_clean1.name), 'Interpreter', 'none');
saveas(fig4, fullfile(fig_dir, 'fig04_clean_window_noisy_psd.png'));

% ---------- Figure 5: delta b and delta kappa ----------
fig5 = figure('Color', 'w');
tl5 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot_delta_shift_map(all_runs, 'fit_img', gca);
title('Image diameter | color = sinking law | shape = frag level');

nexttile;
plot_delta_shift_map(all_runs, 'fit_vol', gca);
title('Volume diameter | color = sinking law | shape = frag level');

title(tl5, 'Relative shift from each law''s no-frag baseline');
saveas(fig5, fullfile(fig_dir, 'fig05_delta_b_delta_kappa_relative.png'));

% ---------- Figure 6: 1 L vs 10 L noise comparison ----------
sample_compare = build_sampled_delta_summary( ...
    all_runs, opts.sample_compare_L, opts.sample_repeats, fit_range_cm, opts.random_seed);

fig6 = figure('Color', 'w');
tl6 = tiledlayout(1, numel(sample_compare), 'TileSpacing', 'compact', 'Padding', 'compact');

for iv = 1:numel(sample_compare)
    nexttile;
    plot_sampled_delta_map(sample_compare(iv), gca);
    title(sprintf('Image diameter sampled | %.0f L | color = law | shape = frag level', ...
        sample_compare(iv).volume_L));
end

title(tl6, sprintf('Sampled delta-space comparison | median of %d repeats', opts.sample_repeats));
saveas(fig6, fullfile(fig_dir, 'fig06_noise_1L_vs_10L_delta_space.png'));

% ---------- Note ----------
note_file = fullfile(note_dir, 'first_pass_figure_note.md');
write_note(note_file, all_runs, law_idx, frag_idx, opts, sample_compare);

out = struct();
out.results = res;
out.note_file = note_file;
out.figure_dir = fig_dir;

fprintf('\nSaved first-pass paper figures to: %s\n', fig_dir);
fprintf('Saved figure note to: %s\n', note_file);

end

function plot_kappa_vs_level(all_runs, field_name)
law_cols = lines(numel(all_runs));
xv = 1:4;
xticks(gca, xv);
xticklabels(gca, {'no-frag', 'weak', 'medium', 'strong'});
xtickangle(gca, 20);

for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    kappav = arrayfun(@(c) c.(field_name).kappa, cases);
    plot(xv, kappav, '-o', 'Color', law_cols(ilaw,:), 'LineWidth', 1.5, ...
        'MarkerSize', 6, 'MarkerFaceColor', law_cols(ilaw,:), ...
        'DisplayName', all_runs(ilaw).law);
end
end

function plot_delta_shift_map(all_runs, field_name, ax)
if nargin < 3 || isempty(ax)
    figure('Color', 'w');
    ax = gca;
end

hold(ax, 'on');
law_cols = lines(numel(all_runs));
markers = {'o', '^', 's', 'd'};

plot(ax, 0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2, ...
    'DisplayName', 'origin');

for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    fit0 = cases(1).(field_name);
    db = zeros(1, numel(cases));
    dk = zeros(1, numel(cases));
    for ic = 1:numel(cases)
        fitc = cases(ic).(field_name);
        db(ic) = fitc.b - fit0.b;
        dk(ic) = fitc.kappa - fit0.kappa;
    end

    plot(ax, db, dk, '-', 'Color', law_cols(ilaw,:), 'LineWidth', 1.3, ...
        'DisplayName', all_runs(ilaw).law);

    for ic = 1:numel(cases)
        plot(ax, db(ic), dk(ic), markers{ic}, ...
            'Color', law_cols(ilaw,:), ...
            'MarkerFaceColor', law_cols(ilaw,:), ...
            'MarkerSize', 7, ...
            'LineStyle', 'none', ...
            'HandleVisibility', 'off');
    end
end

grid(ax, 'on');
xlabel(ax, 'delta b');
ylabel(ax, 'delta kappa');
legend(ax, 'Location', 'best');
text(ax, 0.02, 0.02, ...
    'marker shape: circle = no-frag, triangle = weak, square = medium, diamond = strong', ...
    'Units', 'normalized', 'FontSize', 8, 'Color', [0.25 0.25 0.25], ...
    'Interpreter', 'none');
end

function summary = build_sampled_delta_summary(all_runs, volumes_L, nrep, fit_range_cm, seed0)
summary = struct([]);

for iv = 1:numel(volumes_L)
    summary(iv).volume_L = volumes_L(iv);
    summary(iv).laws = struct([]);

    for ilaw = 1:numel(all_runs)
        cases = all_runs(ilaw).cases;
        fit0 = cases(1).fit_img;
        db = nan(1, numel(cases));
        dk = nan(1, numel(cases));
        good = zeros(1, numel(cases));

        for ic = 1:numel(cases)
            bvals = nan(nrep, 1);
            kvals = nan(nrep, 1);
            D = cases(ic).res.output_data.diam_i(:);
            N = cases(ic).res.output_data.nspec_i(end, :)';

            for ir = 1:nrep
                rng(seed0 + 1000 * iv + 100 * ilaw + 10 * ic + ir);
                fit_s = sample_fit_local(D, N, fit_range_cm, volumes_L(iv));
                if fit_s.ok
                    bvals(ir) = fit_s.b;
                    kvals(ir) = fit_s.kappa;
                end
            end

            db(ic) = median(bvals, 'omitnan') - fit0.b;
            dk(ic) = median(kvals, 'omitnan') - fit0.kappa;
            good(ic) = sum(isfinite(bvals));
        end

        summary(iv).laws(ilaw).law = all_runs(ilaw).law;
        summary(iv).laws(ilaw).db = db;
        summary(iv).laws(ilaw).dk = dk;
        summary(iv).laws(ilaw).good = good;
    end
end
end

function plot_sampled_delta_map(summary, ax)
if nargin < 2 || isempty(ax)
    figure('Color', 'w');
    ax = gca;
end

hold(ax, 'on');
law_cols = lines(numel(summary.laws));
markers = {'o', '^', 's', 'd'};

plot(ax, 0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2, ...
    'DisplayName', 'clean no-frag baseline');

for ilaw = 1:numel(summary.laws)
    db = summary.laws(ilaw).db;
    dk = summary.laws(ilaw).dk;

    plot(ax, db, dk, '-', 'Color', law_cols(ilaw,:), 'LineWidth', 1.3, ...
        'DisplayName', summary.laws(ilaw).law);

    for ic = 1:numel(db)
        plot(ax, db(ic), dk(ic), markers{ic}, ...
            'Color', law_cols(ilaw,:), ...
            'MarkerFaceColor', law_cols(ilaw,:), ...
            'MarkerSize', 7, ...
            'LineStyle', 'none', ...
            'HandleVisibility', 'off');
    end
end

grid(ax, 'on');
xlabel(ax, 'delta b sampled');
ylabel(ax, 'delta kappa sampled');
legend(ax, 'Location', 'best');
text(ax, 0.02, 0.02, ...
    'marker shape: circle = no-frag, triangle = weak, square = medium, diamond = strong', ...
    'Units', 'normalized', 'FontSize', 8, 'Color', [0.25 0.25 0.25], ...
    'Interpreter', 'none');
end

function plot_clean_pair(fit0, fit1, show_legend)
loglog(fit0.Dall, fit0.Nall, '-', 'Color', fit0.color, 'LineWidth', 1.4, ...
    'DisplayName', fit0.level_name);
loglog(fit1.Dall, fit1.Nall, '-', 'Color', fit1.color, 'LineWidth', 1.4, ...
    'DisplayName', fit1.level_name);
loglog(fit0.Dfit, fit0.Nfit, '--', 'Color', fit0.color, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
loglog(fit1.Dfit, fit1.Nfit, '--', 'Color', fit1.color, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Image diameter (cm)');
ylabel('N(D)');
if show_legend
    legend('Location', 'best');
end
end

function plot_window_pair(fit0, fit1, fit_range_cm)
plot(fit0.x, fit0.resid, '-o', 'Color', fit0.color, 'LineWidth', 1.2, ...
    'MarkerSize', 4, 'DisplayName', sprintf('%s | kappa=%.3g', fit0.level_name, fit0.kappa));
plot(fit1.x, fit1.resid, '-o', 'Color', fit1.color, 'LineWidth', 1.2, ...
    'MarkerSize', 4, 'DisplayName', sprintf('%s | kappa=%.3g', fit1.level_name, fit1.kappa));
plot(fit0.x, fit0.resid_quad, '--', 'Color', fit0.color, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
plot(fit1.x, fit1.resid_quad, '--', 'Color', fit1.color, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on;
xlabel('log10(D)');
ylabel('Residual');
title(sprintf('Fit window %.3g to %.3g cm', fit_range_cm(1), fit_range_cm(2)));
legend('Location', 'best');
end

function plot_noisy_pair(fit0, fit0s, fit1, fit1s, sample_L)
loglog(fit0.Dall, fit0.Nall, '-', 'Color', fit0.color, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
loglog(fit1.Dall, fit1.Nall, '-', 'Color', fit1.color, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');

loglog(fit0s.Dall, fit0s.Nall, 'o', 'Color', fit0.color, 'MarkerSize', 4, ...
    'DisplayName', sprintf('%s sampled', fit0.level_name));
loglog(fit1s.Dall, fit1s.Nall, 'o', 'Color', fit1.color, 'MarkerSize', 4, ...
    'DisplayName', sprintf('%s sampled', fit1.level_name));

if fit0s.ok
    loglog(fit0s.Dfit, fit0s.Nfit, '--', 'Color', fit0.color, 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end
if fit1s.ok
    loglog(fit1s.Dfit, fit1s.Nfit, '--', 'Color', fit1.color, 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end

grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Image diameter (cm)');
ylabel('N(D) sampled');
title(sprintf('Poisson sampled | V=%.2g L', sample_L));
legend('Location', 'best');
end

function [law_idx, frag_idx] = strongest_case_idx(all_runs)
law_idx = 1;
frag_idx = 2;
best_score = -Inf;

for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    k0 = cases(1).fit_img.kappa;
    for ic = 2:numel(cases)
        score = abs(cases(ic).fit_img.kappa - k0);
        if score > best_score
            best_score = score;
            law_idx = ilaw;
            frag_idx = ic;
        end
    end
end
end

function idx = pick_reference_law(all_runs)
idx = 1;
for i = 1:numel(all_runs)
    if strcmpi(all_runs(i).law, 'kriest_8')
        idx = i;
        return
    end
end
end

function fit_s = sample_fit_local(D, N, fit_rng, sample_L)
[~, dD] = bin_edges_from_centers_local(D);
N_L = max(N(:), 0) * 1e3;
lambda = N_L .* dD .* sample_L;
counts = poisson_sample_local(lambda);
N_obs = counts ./ max(sample_L * dD, eps);
fit_s = compute_psd_diagnostics(D, N_obs, fit_rng);
end

function [edges, dD] = bin_edges_from_centers_local(D)
D = D(:);
n = numel(D);
edges = zeros(n + 1, 1);

if n < 2
    edges(1) = D(1) / 1.2;
    edges(2) = D(1) * 1.2;
    dD = diff(edges);
    return
end

edges(2:n) = sqrt(D(1:n-1) .* D(2:n));
edges(1) = D(1)^2 / edges(2);
edges(n + 1) = D(n)^2 / edges(n);
dD = diff(edges);
dD(dD <= 0) = eps;
end

function k = poisson_sample_local(lambda)
lambda = lambda(:);
k = zeros(size(lambda));
for i = 1:numel(lambda)
    L = max(lambda(i), 0);
    if L < 30
        p = 1.0;
        n = 0;
        lim = exp(-L);
        while p > lim
            n = n + 1;
            p = p * rand();
        end
        k(i) = n - 1;
    else
        k(i) = max(0, round(L + sqrt(L) * randn()));
    end
end
end

function fits = gather_fits(all_runs, field_name)
fits = struct([]);
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        fit = cases(ic).(field_name);
        if ~isempty(fit) && isfield(fit, 'ok') && fit.ok
            if isempty(fits)
                fits = fit;
            else
                fits(end + 1) = fit; %#ok<AGROW>
            end
        end
    end
end
end

function write_note(note_file, all_runs, law_idx, frag_idx, opts, sample_compare)
fid = fopen(note_file, 'w');
if fid < 0
    return
end

law_name = all_runs(law_idx).law;
case_name = all_runs(law_idx).cases(frag_idx).name;
fit0 = all_runs(law_idx).cases(1).fit_img;
fit1 = all_runs(law_idx).cases(frag_idx).fit_img;
dk = fit1.kappa - fit0.kappa;

fprintf(fid, '# First-pass Figure Note\n\n');
fprintf(fid, 'The strongest figure for the paper is currently `fig05_delta_b_delta_kappa_relative.png`.\n\n');
fprintf(fid, 'Why it looks strongest:\n\n');
fprintf(fid, '- it removes the large baseline offsets among sinking laws\n');
fprintf(fid, '- it shows the fragmentation response as a shift from each law''s own no-frag reference\n');
fprintf(fid, '- it makes the image-space contrast easier to read than the raw `(b, kappa)` plot\n\n');
fprintf(fid, 'The strongest single image-space contrast in this run is `%s | %s` relative to its no-frag baseline.\n\n', ...
    law_name, case_name);
fprintf(fid, '- baseline image-space `kappa` = `%.4g`\n', fit0.kappa);
fprintf(fid, '- case image-space `kappa` = `%.4g`\n', fit1.kappa);
fprintf(fid, '- shift in `kappa` = `%.4g`\n\n', dk);
fprintf(fid, 'Support figures:\n\n');
fprintf(fid, '- `fig01_b_kappa_separability.png` is still useful as the raw metric-space view.\n');
fprintf(fid, '- `fig02_baseline_psd_fit_residuals.png` is good for showing that slope-only looks too simple.\n');
fprintf(fid, '- `fig03_kappa_vs_fragmentation_strength.png` is good for trend summary.\n');
fprintf(fid, '- `fig04_clean_window_noisy_psd.png` is good for showing what may be lost once the PSD is windowed and sampled.\n');
fprintf(fid, '- `fig05_delta_b_delta_kappa_relative.png` is good for removing law-to-law baseline offsets.\n');
fprintf(fid, '- `fig06_noise_1L_vs_10L_delta_space.png` is good for showing that larger sample volume stabilizes the sampled metric space.\n\n');
fprintf(fid, 'Note:\n\n');
fprintf(fid, '- this note is written by the paper-figure script\n');
fprintf(fid, '- UVP-like sampling was `%s`\n', mat2str(logical(opts.apply_uvp)));
if isfield(opts, 'uvp_sample_L')
    fprintf(fid, '- sample volume = `%.4g L`\n', opts.uvp_sample_L);
end
if nargin >= 6 && ~isempty(sample_compare)
    fprintf(fid, '- extra sampled delta figure uses `%d` repeats at `%.0f` and `%.0f` L\n', ...
        opts.sample_repeats, sample_compare(1).volume_L, sample_compare(end).volume_L);
end

fclose(fid);
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
