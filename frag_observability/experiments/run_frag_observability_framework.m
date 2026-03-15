function out = run_frag_observability_framework(opts)
% run_frag_observability_framework
% Small paper experiment set for frag observability.

if nargin < 1 || isempty(opts)
    opts = struct();
end

opts = set_default(opts, 'fit_range_cm', [0.01, 0.6]);
opts = set_default(opts, 'plot_range_cm', [0.005, 2.0]);
opts = set_default(opts, 'sinking_laws', {'current', 'kriest_8', 'kriest_9', 'siegel_2025'});
opts = set_default(opts, 'scale_scan', [0.8, 1.0, 1.2]);
opts = set_default(opts, 'frag_eps', [1e-8, 1e-6, 1e-4]);
opts = set_default(opts, 'frag_names', {'weak', 'medium', 'strong'});
opts = set_default(opts, 'apply_uvp', false);
opts = set_default(opts, 'uvp_sample_L', 1.0);
opts = set_default(opts, 'random_seed', 1);

exp_dir = fileparts(mfilename('fullpath'));
frag_root = fileparts(exp_dir);
repo_root = fileparts(frag_root);
model_root = fullfile(repo_root, 'coag_model_final');
diag_root = fullfile(frag_root, 'diagnostics');
fig_dir = fullfile(frag_root, 'figures');
tbl_dir = fullfile(frag_root, 'tables');

if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(tbl_dir, 'dir'), mkdir(tbl_dir); end

run(fullfile(model_root, 'setup_paths.m'));
addpath(diag_root);
rehash;

rng(opts.random_seed);

base = make_base_cfg();
all_runs = struct([]);
rows = struct([]);
base_rows = struct([]);

case_cols = [ ...
    0.10 0.10 0.10; ...
    0.16 0.45 0.78; ...
    0.86 0.40 0.12; ...
    0.73 0.16 0.18];
case_markers = {'o', '^', 's', 'd'};

for ilaw = 1:numel(opts.sinking_laws)
    law = opts.sinking_laws{ilaw};
    fprintf('\n=== law: %s ===\n', law);

    base_law = base.copy();
    base_law.sinking_law = law;

    scan = struct([]);
    for is = 1:numel(opts.scale_scan)
        cfg = base_law.copy();
        cfg.sinking_scale = opts.scale_scan(is);
        cfg.enable_disagg = false;
        cfg.disagg_mode = 'legacy';

        sim = CoagulationSimulation(cfg);
        res = sim.run();

        fit_img = compute_psd_diagnostics( ...
            res.output_data.diam_i(:), ...
            res.output_data.nspec_i(end, :)', ...
            opts.fit_range_cm);

        scan(is).cfg = cfg;
        scan(is).res = res;
        scan(is).fit_img = fit_img;
        scan(is).scale = opts.scale_scan(is);

        fprintf('scale=%.2f | image R2=%.4f | b=%.4f\n', ...
            opts.scale_scan(is), fit_img.r2, fit_img.b);
    end

    [~, best_i] = max(arrayfun(@(s) s.fit_img.r2, scan));
    best_scale = scan(best_i).scale;

    base_rows(ilaw).sinking_law = law;
    base_rows(ilaw).best_scale = best_scale;
    base_rows(ilaw).r2_no_frag_image = scan(best_i).fit_img.r2;
    base_rows(ilaw).b_no_frag_image = scan(best_i).fit_img.b;

    fprintf('picked scale %.2f for %s\n', best_scale, law);

    cases = struct([]);
    level_names = [{'no_frag'}, opts.frag_names];
    level_eps = [NaN, opts.frag_eps];

    for ic = 1:numel(level_names)
        cfg = base_law.copy();
        cfg.sinking_scale = best_scale;

        if ic == 1
            cfg.enable_disagg = false;
            cfg.disagg_mode = 'legacy';
        else
            cfg.enable_disagg = true;
            cfg.disagg_mode = 'operator_split';
            cfg.disagg_epsilon = opts.frag_eps(ic - 1);
            cfg.disagg_outer_dt = 1/24;
            cfg.disagg_dmax_cm = [];
        end

        sim = CoagulationSimulation(cfg);
        res = sim.run();

        [dmax_vol_cm, dmax_img_cm] = get_dmax_markers(cfg, res);

        fit_img = compute_psd_diagnostics( ...
            res.output_data.diam_i(:), ...
            res.output_data.nspec_i(end, :)', ...
            opts.fit_range_cm);
        fit_img = add_fit_meta(fit_img, law, 'image', level_names{ic}, ...
            case_cols(ic,:), case_markers{ic});

        fit_vol = compute_psd_diagnostics( ...
            res.output_data.diam_v(:), ...
            res.output_data.nspec_v(end, :)', ...
            opts.fit_range_cm);
        fit_vol = add_fit_meta(fit_vol, law, 'volume', level_names{ic}, ...
            case_cols(ic,:), case_markers{ic});

        sampled_img = struct('ok', false);
        sampled_vol = struct('ok', false);
        if opts.apply_uvp
            sampled_img = sample_fit( ...
                res.output_data.diam_i(:), ...
                res.output_data.nspec_i(end, :)', ...
                opts.fit_range_cm, opts.uvp_sample_L);
            sampled_img = add_fit_meta(sampled_img, law, 'image_sampled', ...
                level_names{ic}, case_cols(ic,:), case_markers{ic});

            sampled_vol = sample_fit( ...
                res.output_data.diam_v(:), ...
                res.output_data.nspec_v(end, :)', ...
                opts.fit_range_cm, opts.uvp_sample_L);
            sampled_vol = add_fit_meta(sampled_vol, law, 'volume_sampled', ...
                level_names{ic}, case_cols(ic,:), case_markers{ic});
        end

        cases(ic).name = level_names{ic};
        cases(ic).eps = level_eps(ic);
        cases(ic).cfg = cfg;
        cases(ic).res = res;
        cases(ic).fit_img = fit_img;
        cases(ic).fit_vol = fit_vol;
        cases(ic).fit_img_sampled = sampled_img;
        cases(ic).fit_vol_sampled = sampled_vol;
        cases(ic).dmax_vol_cm = dmax_vol_cm;
        cases(ic).dmax_img_cm = dmax_img_cm;

        rows = append_row(rows, law, best_scale, level_names{ic}, level_eps(ic), ...
            'image', false, NaN, fit_img, dmax_vol_cm, dmax_img_cm);
        rows = append_row(rows, law, best_scale, level_names{ic}, level_eps(ic), ...
            'volume', false, NaN, fit_vol, dmax_vol_cm, dmax_img_cm);

        if opts.apply_uvp
            rows = append_row(rows, law, best_scale, level_names{ic}, level_eps(ic), ...
                'image', true, opts.uvp_sample_L, sampled_img, dmax_vol_cm, dmax_img_cm);
            rows = append_row(rows, law, best_scale, level_names{ic}, level_eps(ic), ...
                'volume', true, opts.uvp_sample_L, sampled_vol, dmax_vol_cm, dmax_img_cm);
        end

        fprintf('%s | eps=%g | b_img=%.4f | kappa_img=%.4g | b_vol=%.4f | kappa_vol=%.4g\n', ...
            level_names{ic}, level_eps(ic), fit_img.b, fit_img.kappa, fit_vol.b, fit_vol.kappa);
    end

    all_runs(ilaw).law = law;
    all_runs(ilaw).best_scale = best_scale;
    all_runs(ilaw).cases = cases;
end

summary_tbl = struct2table(rows);
base_tbl = struct2table(base_rows);

writetable(summary_tbl, fullfile(tbl_dir, 'frag_observability_case_summary.csv'));
writetable(base_tbl, fullfile(tbl_dir, 'frag_observability_baseline_scales.csv'));
save(fullfile(tbl_dir, 'frag_observability_case_summary.mat'), ...
    'all_runs', 'summary_tbl', 'base_tbl', 'opts');

make_figures(all_runs, fig_dir, opts);

out = struct();
out.all_runs = all_runs;
out.summary_tbl = summary_tbl;
out.base_tbl = base_tbl;
out.opts = opts;

fprintf('\nSaved tables to: %s\n', tbl_dir);
fprintf('Saved figures to: %s\n', fig_dir);

end

function cfg = make_base_cfg()
cfg = SimulationConfig();
cfg.t_final = 200;
cfg.delta_t = 1.0;
cfg.n_sections = 35;
cfg.enable_pp = true;
cfg.pp_bin = 1;
cfg.pp_source = 1e-8;
cfg.enable_coag = true;
cfg.enable_sinking = true;
cfg.enable_linear = true;
cfg.growth = 0;
cfg.box_depth = 2000;
cfg.r_to_rg = 1.6;
cfg.sinking_size = 'image';
end

function fit = add_fit_meta(fit, law, space, level_name, color, marker)
fit.name = sprintf('%s | %s | %s', law, space, level_name);
fit.group = law;
fit.color = color;
fit.marker = marker;
fit.space = space;
fit.level_name = level_name;
end

function rows = append_row(rows, law, best_scale, level_name, eps_frag, space, sampled, sample_L, fit, dmax_vol_cm, dmax_img_cm)
row = struct();
row.sinking_law = string(law);
row.best_scale = best_scale;
row.frag_level = string(level_name);
row.frag_eps = eps_frag;
row.diameter_space = string(space);
row.sampled = logical(sampled);
row.sample_volume_L = sample_L;
row.b = fit.b;
row.kappa = fit.kappa;
row.r2 = fit.r2;
row.alpha = fit.alpha;
row.c0 = fit.c0;
row.c1 = fit.c1;
row.c2 = fit.c2;
row.n_fit = fit.n_fit;
row.fit_dmin_cm = fit.fit_range(1);
row.fit_dmax_cm = fit.fit_range(2);
row.dmax_vol_cm = dmax_vol_cm;
row.dmax_img_cm = dmax_img_cm;

if isempty(rows)
    rows = row;
else
    rows(end + 1) = row; %#ok<AGROW>
end
end

function [dmax_vol_cm, dmax_img_cm] = get_dmax_markers(cfg, res)
dmax_vol_cm = NaN;
dmax_img_cm = NaN;

if ~isprop(cfg, 'enable_disagg') || ~cfg.enable_disagg
    return
end
if ~isprop(cfg, 'disagg_mode') || ~strcmpi(string(cfg.disagg_mode), 'operator_split')
    return
end

if isprop(cfg, 'disagg_dmax_cm') && ~isempty(cfg.disagg_dmax_cm)
    dmax_vol_cm = cfg.disagg_dmax_cm;
elseif isprop(cfg, 'disagg_epsilon') && ~isempty(cfg.disagg_epsilon)
    dmax_vol_cm = 0.1 * cfg.disagg_C * cfg.disagg_epsilon^(-cfg.disagg_gamma);
end

if ~isfinite(dmax_vol_cm)
    return
end

Dv = res.output_data.diam_v(:);
Di = res.output_data.diam_i(:);
if ~isempty(Dv) && ~isempty(Di) && numel(Dv) == numel(Di)
    [~, idx] = min(abs(Dv - dmax_vol_cm));
    dmax_img_cm = Di(idx);
end
end

function sampled_fit = sample_fit(D, N, fit_rng, sample_L)
[~, dD] = bin_edges_from_centers(D);
N_L = max(N(:), 0) * 1e3;
lambda = N_L .* dD .* sample_L;
counts = poisson_sample(lambda);
N_obs = counts ./ max(sample_L * dD, eps);
sampled_fit = compute_psd_diagnostics(D, N_obs, fit_rng);
sampled_fit.sample_counts = counts;
sampled_fit.sample_volume_L = sample_L;
end

function [edges, dD] = bin_edges_from_centers(D)
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

function k = poisson_sample(lambda)
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

function make_figures(all_runs, fig_dir, opts)

fig1 = figure('Color', 'w');
tl1 = tiledlayout(numel(all_runs), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ilaw = 1:numel(all_runs)
    nexttile;
    fits_img = [all_runs(ilaw).cases.fit_img];
    plot_psd_with_fit(fits_img, gca);
    xlim(opts.plot_range_cm);
    title(sprintf('%s | image', all_runs(ilaw).law), 'Interpreter', 'none');

    nexttile;
    fits_vol = [all_runs(ilaw).cases.fit_vol];
    plot_psd_with_fit(fits_vol, gca);
    xlim(opts.plot_range_cm);
    title(sprintf('%s | volume', all_runs(ilaw).law), 'Interpreter', 'none');
end
title(tl1, 'Frag observability PSD fits');
saveas(fig1, fullfile(fig_dir, 'frag_observability_psd_with_fit.png'));

fig2 = figure('Color', 'w');
tl2 = tiledlayout(numel(all_runs), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ilaw = 1:numel(all_runs)
    nexttile;
    fits_img = [all_runs(ilaw).cases.fit_img];
    plot_psd_residuals(fits_img, gca);
    title(sprintf('%s | image', all_runs(ilaw).law), 'Interpreter', 'none');

    nexttile;
    fits_vol = [all_runs(ilaw).cases.fit_vol];
    plot_psd_residuals(fits_vol, gca);
    title(sprintf('%s | volume', all_runs(ilaw).law), 'Interpreter', 'none');
end
title(tl2, 'Frag observability residuals');
saveas(fig2, fullfile(fig_dir, 'frag_observability_residuals.png'));

fig3 = figure('Color', 'w');
tl3 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot_b_kappa_separability(gather_fits(all_runs, 'fit_img'), gca, false);
title('image');
nexttile;
plot_b_kappa_separability(gather_fits(all_runs, 'fit_vol'), gca, false);
title('volume');
title(tl3, '(b, kappa) separability');
saveas(fig3, fullfile(fig_dir, 'frag_observability_b_kappa.png'));

if opts.apply_uvp
    fig4 = figure('Color', 'w');
    tl4 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    plot_b_kappa_separability(gather_fits(all_runs, 'fit_img_sampled'), gca, false);
    title(sprintf('image sampled | V=%.2g L', opts.uvp_sample_L));
    nexttile;
    plot_b_kappa_separability(gather_fits(all_runs, 'fit_vol_sampled'), gca, false);
    title(sprintf('volume sampled | V=%.2g L', opts.uvp_sample_L));
    title(tl4, '(b, kappa) separability after Poisson sample');
    saveas(fig4, fullfile(fig_dir, 'frag_observability_b_kappa_sampled.png'));
end

end

function fits = gather_fits(all_runs, field_name)
fits = struct([]);
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        if isfield(cases(ic), field_name)
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
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
