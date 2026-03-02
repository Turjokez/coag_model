% run_powerlaw_frag_levels.m
% Main checks:
% 1) get a near power-law baseline
% 2) run no-frag steady case
% 3) add frag at 1e-8, 1e-6, 1e-4
% Optional: add one extra strong test for sensitivity.

clear; close all; clc;
setup_paths

% ---------- output folder ----------
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
fig_dir = fullfile(proj_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% ---------- experiment knobs ----------
fit_range_cm = [0.01, 0.6];           % fit window for log-log slope
plot_range_cm = [0.005, 2.0];         % show broad image-diameter range
eps_list_core = [1e-8, 1e-6, 1e-4];   % low / medium / high
eps_list_extra = [];                   % optional, e.g. [1e-2]
eps_list = [eps_list_core, eps_list_extra];
scale_list = [0.8, 1.0, 1.2];         % quick search for near power-law
sinking_law = 'kriest_8';
steady_t0 = 20;                       % ignore startup when comparing slopes
uvp_sample_L = 1.0;                   % synthetic UVP sample volume [L]

% ---------- baseline config ----------
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

base.sinking_law = sinking_law;
base.sinking_size = 'image';

% ---------- Step 1: pick best no-frag power-law ----------
fprintf('\n=== Step 1: power-law search (no frag) ===\n');
best = struct('r2', -Inf);
scan = struct([]);

for i = 1:numel(scale_list)
    cfg = base.copy();
    cfg.sinking_scale = scale_list(i);
    cfg.enable_disagg = false;
    cfg.disagg_mode = 'legacy';

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    D = res.output_data.diam_i(:);
    N = res.output_data.nspec_i(end, :)';
    fit = fit_powerlaw_loglog(D, N, fit_range_cm);

    scan(i).scale = scale_list(i);
    scan(i).fit = fit;
    scan(i).res = res;

    fprintf('scale=%.2f | R2=%.4f | slope=%.4f | b=%.4f\n', ...
        scale_list(i), fit.r2, fit.slope, fit.b);

    if fit.r2 > best.r2
        best = fit;
        best_idx = i; %#ok<NASGU>
    end
end

[~, best_i] = max(arrayfun(@(s) s.fit.r2, scan));
best_scale = scan(best_i).scale;
res_no_frag = scan(best_i).res;
fit_no_frag = scan(best_i).fit;

fprintf('\nSelected no-frag baseline: law=%s scale=%.2f (R2=%.4f)\n', ...
    sinking_law, best_scale, fit_no_frag.r2);

% ---------- Step 2/3: frag levels ----------
fprintf('\n=== Step 2/3: add fragmentation by epsilon ===\n');
fprintf('Core eps: %s | extra stress eps: %s\n', ...
    mat2str(eps_list_core), mat2str(eps_list_extra));
cases = struct([]);

% case 1: no frag
cases(1).name = 'no_frag';
cases(1).eps = NaN;
cases(1).cfg = base.copy();
cases(1).cfg.sinking_scale = best_scale;
cases(1).cfg.enable_disagg = false;
cases(1).cfg.disagg_mode = 'legacy';
cases(1).res = res_no_frag;
cases(1).fit = fit_no_frag;
cases(1).dmax_cm = NaN;
cases(1).dmax_cm_vol = NaN;
cases(1).dmax_cm_img = NaN;
cases(1).is_core = true;

% cases 2-4: operator-split frag
for j = 1:numel(eps_list)
    idx = j + 1;
    cfg = base.copy();
    cfg.sinking_scale = best_scale;
    cfg.enable_disagg = true;
    cfg.disagg_mode = 'operator_split';
    cfg.disagg_epsilon = eps_list(j);
    cfg.disagg_outer_dt = 1/24;
    cfg.disagg_dmax_cm = [];

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    D = res.output_data.diam_i(:);
    N = res.output_data.nspec_i(end, :)';
    fit = fit_powerlaw_loglog(D, N, fit_range_cm);

    dmax_cm_vol = 0.1 * cfg.disagg_C * cfg.disagg_epsilon^(-cfg.disagg_gamma);
    Dv = res.output_data.diam_v(:);
    Di = res.output_data.diam_i(:);
    dmax_cm_img = NaN;
    if ~isempty(Dv) && ~isempty(Di) && numel(Dv) == numel(Di)
        [~, i_dmax] = min(abs(Dv - dmax_cm_vol));
        dmax_cm_img = Di(i_dmax);
    end

    is_core = ismember(eps_list(j), eps_list_core);
    if is_core
        cases(idx).name = sprintf('frag_eps_%0.0e', eps_list(j));
    else
        cases(idx).name = sprintf('frag_eps_%0.0e_stress', eps_list(j));
    end
    cases(idx).eps = eps_list(j);
    cases(idx).cfg = cfg;
    cases(idx).res = res;
    cases(idx).fit = fit;
    cases(idx).dmax_cm = dmax_cm_vol;
    cases(idx).dmax_cm_vol = dmax_cm_vol;
    cases(idx).dmax_cm_img = dmax_cm_img;
    cases(idx).is_core = is_core;

    fprintf('eps=%0.0e | Dmax(vol)=%.3f cm | Dmax(img)=%.3f cm | R2=%.4f | slope=%.4f | b=%.4f\n', ...
        eps_list(j), dmax_cm_vol, dmax_cm_img, fit.r2, fit.slope, fit.b);
end

% ---------- fit slope over time ----------
for ic = 1:numel(cases)
    t = cases(ic).res.time(:);
    D = cases(ic).res.output_data.diam_i(:);
    Ns = cases(ic).res.output_data.nspec_i;
    b_t = nan(numel(t),1);
    r2_t = nan(numel(t),1);
    for k = 1:numel(t)
        fit_k = fit_powerlaw_loglog(D, Ns(k,:)', fit_range_cm);
        b_t(k) = fit_k.b;
        r2_t(k) = fit_k.r2;
    end
    cases(ic).b_t = b_t;
    cases(ic).r2_t = r2_t;

    Nend = cases(ic).res.output_data.nspec_i(end, :)';
    if numel(Nend) >= 2
        cases(ic).tail_ratio_last = Nend(end) / max(Nend(end-1), eps);
    else
        cases(ic).tail_ratio_last = NaN;
    end
end

% ---------- Figure 1: final spectra + fits ----------
fig1 = figure('Color','w'); hold on;
cols = lines(numel(cases));
for ic = 1:numel(cases)
    D = cases(ic).res.output_data.diam_i(:);
    N = cases(ic).res.output_data.nspec_i(end, :)';
    fit = cases(ic).fit;
    label = case_label(cases(ic));
    loglog(D, N, '-', 'LineWidth', 1.4, 'Color', cols(ic,:), 'DisplayName', label);
    if fit.ok
        Df = fit.Dfit;
        Nf = fit.a * (Df.^(-fit.b));
        loglog(Df, Nf, '--', 'LineWidth', 1.1, 'Color', cols(ic,:), 'HandleVisibility','off');
    end
end
grid on;
xlabel('Image diameter (cm)');
ylabel('Number spectrum N');
title('Final size spectrum: no-frag vs fragmentation levels');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(plot_range_cm);
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'powerlaw_frag_levels_final_spectra.png'));

% ---------- Figure 2: slope b over time (steady window) ----------
fig2 = figure('Color','w'); hold on;
for ic = 1:numel(cases)
    t = cases(ic).res.time(:);
    mask = t >= steady_t0;
    plot(t(mask), cases(ic).b_t(mask), 'LineWidth', 1.4, 'Color', cols(ic,:), ...
        'DisplayName', case_label(cases(ic)));
end
grid on;
xlabel('Time (days)');
ylabel('Power-law exponent b');
title(sprintf('Power-law exponent over time (t >= %g d)', steady_t0));
legend('Location','best');
saveas(fig2, fullfile(fig_dir, 'powerlaw_frag_levels_b_timeseries_steady.png'));

% ---------- Figure 3: final residuals ----------
fig3 = figure('Color','w'); hold on;
for ic = 1:numel(cases)
    fit = cases(ic).fit;
    if ~fit.ok
        continue
    end
    plot(log10(fit.Dsel), fit.resid, '-o', 'LineWidth', 1.1, 'MarkerSize', 3, ...
        'Color', cols(ic,:), 'DisplayName', case_label(cases(ic)));
end
yline(0, 'k--', 'LineWidth', 1.0);
grid on;
xlabel('log10(D)');
ylabel('Residual: log10(N) - fit');
title('Final residual pattern in fit window');
legend('Location','best');
saveas(fig3, fullfile(fig_dir, 'powerlaw_frag_levels_residuals.png'));

% ---------- Figure 4: frag/no-frag ratio ----------
fig4 = figure('Color','w'); hold on;
D0 = cases(1).res.output_data.diam_i(:);
N0 = cases(1).res.output_data.nspec_i(end, :)';
for ic = 2:numel(cases)
    D = cases(ic).res.output_data.diam_i(:);
    N = cases(ic).res.output_data.nspec_i(end, :)';
    mask = (N0 > 1e-3 * max(N0)) & (D >= plot_range_cm(1)) & (D <= plot_range_cm(2));
    ratio = nan(size(N));
    ratio(mask) = N(mask) ./ N0(mask);
    loglog(D, ratio, '-', 'LineWidth', 1.4, 'Color', cols(ic,:), ...
        'DisplayName', case_label(cases(ic)));
end
yline(1.0, 'k--', 'LineWidth', 1.0, 'DisplayName','no frag ref');
grid on;
xlabel('Image diameter (cm)');
ylabel('N_{frag} / N_{no-frag}');
title('How fragmentation changes shape (final time)');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(plot_range_cm);
legend('Location','best');
saveas(fig4, fullfile(fig_dir, 'powerlaw_frag_levels_ratio_vs_nofrag.png'));

% ---------- Figure 5: UVP-like sampled points (Poisson) ----------
% Use no-frag and one frag case for a direct visual check.
fig5 = figure('Color','w');
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;

ic_eps_1e4 = find(arrayfun(@(c) ~isnan(c.eps) && abs(c.eps - 1e-4) < 1e-14, cases), 1, 'first');
if isempty(ic_eps_1e4)
    frag_ids = find(~isnan([cases.eps]));
    ic_eps_1e4 = frag_ids(end);
end
case_uvp = [1, ic_eps_1e4];
uvp_cols = [0.1 0.1 0.8; 0.8 0.2 0.1];
for iu = 1:numel(case_uvp)
    ic = case_uvp(iu);
    D = cases(ic).res.output_data.diam_i(:);
    N = cases(ic).res.output_data.nspec_i(end, :)';
    [edges, dD] = bin_edges_from_centers(D);

    % model mean counts in each bin (Poisson lambda)
    N_L = N * 1e3; % #/L/cm
    lambda = N_L .* dD .* uvp_sample_L;
    counts = poisson_sample(lambda);

    % sampled spectrum estimate
    N_obs_L = counts ./ max(uvp_sample_L * dD, eps);
    one_particle = 1 ./ max(uvp_sample_L * dD, eps);

    % valid data for plot and fit
    ok = counts > 0 & D >= plot_range_cm(1) & D <= plot_range_cm(2);
    Dp = D(ok);
    Np = N_obs_L(ok);

    if isempty(Dp)
        continue
    end

    fit_uvp = fit_powerlaw_loglog(Dp, Np, fit_range_cm);
    cases(ic).uvp_fit = fit_uvp;
    cases(ic).uvp_D = Dp;
    cases(ic).uvp_N = Np;

    % smooth model + sampled points + one-particle line
    loglog(D, N_L, '-', 'Color', uvp_cols(iu,:), 'LineWidth', 1.2, ...
        'DisplayName', [case_label(cases(ic)) ' (model)']);
    loglog(Dp, Np, 'o', 'Color', uvp_cols(iu,:), 'MarkerSize', 4, ...
        'DisplayName', [case_label(cases(ic)) ' (sampled)']);
    loglog(D, one_particle, '--', 'Color', uvp_cols(iu,:), 'LineWidth', 1.0, ...
        'HandleVisibility','off');

    if fit_uvp.ok
        Df = fit_uvp.Dfit;
        Nf = fit_uvp.a * (Df.^(-fit_uvp.b));
        loglog(Df, Nf, ':', 'Color', uvp_cols(iu,:), 'LineWidth', 1.2, ...
            'HandleVisibility','off');
    end
end
grid on;
xlabel('Image diameter (cm)');
ylabel('N (counts/L/cm)');
title('UVP-like sampled points + one-particle line');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(plot_range_cm);
legend('Location','best');

nexttile; hold on;
for iu = 1:numel(case_uvp)
    ic = case_uvp(iu);
    if ~isfield(cases(ic), 'uvp_fit') || ~cases(ic).uvp_fit.ok
        continue
    end
    x = log10(cases(ic).uvp_fit.Dsel);
    r = cases(ic).uvp_fit.resid;
    plot(x, r, '-o', 'Color', uvp_cols(iu,:), 'LineWidth', 1.1, 'MarkerSize', 4, ...
        'DisplayName', case_label(cases(ic)));
end
yline(0, 'k--', 'LineWidth', 1.0);
grid on;
xlabel('log10(D)');
ylabel('Residual (sampled)');
title('Residual pattern from UVP-like sampled points');
legend('Location','best');

title(tl, sprintf('Synthetic UVP check (V=%.1f L)', uvp_sample_L));
saveas(fig5, fullfile(fig_dir, 'powerlaw_frag_levels_uvp_like_sampled.png'));

% ---------- Figure 6: slope vs epsilon (core frag levels) ----------
core_ids = find(arrayfun(@(c) ~isnan(c.eps) && c.is_core, cases));
if ~isempty(core_ids)
    eps_core = [cases(core_ids).eps];
    [eps_core, ord] = sort(eps_core);
    core_ids = core_ids(ord);

    b_final = arrayfun(@(c) c.fit.b, cases(core_ids));
    b_med = nan(size(core_ids));
    b_min = nan(size(core_ids));
    b_max = nan(size(core_ids));
    dmax_core = arrayfun(@(c) c.dmax_cm_vol, cases(core_ids));
    for ii = 1:numel(core_ids)
        ic = core_ids(ii);
        tt = cases(ic).res.time(:);
        mask = tt >= steady_t0;
        bvals = cases(ic).b_t(mask);
        b_med(ii) = median(bvals, 'omitnan');
        b_min(ii) = min(bvals);
        b_max(ii) = max(bvals);
    end

    fig6 = figure('Color','w');
    yyaxis left; hold on;
    semilogx(eps_core, b_final, '-o', 'LineWidth', 1.4, 'MarkerSize', 6, ...
        'DisplayName', 'b (final)');
    semilogx(eps_core, b_med, '-s', 'LineWidth', 1.4, 'MarkerSize', 6, ...
        'DisplayName', sprintf('b (t >= %g d median)', steady_t0));
    for ii = 1:numel(eps_core)
        plot([eps_core(ii), eps_core(ii)], [b_min(ii), b_max(ii)], '-', ...
            'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'HandleVisibility','off');
    end
    yline(cases(1).fit.b, 'k--', 'LineWidth', 1.0, ...
        'DisplayName', 'no-frag b (final)');
    ylabel('Power-law exponent b');

    yyaxis right;
    semilogx(eps_core, dmax_core, ':d', 'LineWidth', 1.2, 'MarkerSize', 6, ...
        'DisplayName', 'D_{max} (cm)');
    ylabel('D_{max} volume-diam (cm)');

    grid on;
    xlabel('Fragmentation epsilon');
    title('Slope response vs epsilon (core levels)');
    legend('Location', 'best');
    saveas(fig6, fullfile(fig_dir, 'powerlaw_frag_levels_slope_vs_epsilon.png'));
end

% ---------- Figure 7: final spectra + Dmax markers ----------
% The model Dmax formula is in volume-diameter space. For this image-diameter
% plot, map Dmax to the nearest image-diameter bin before drawing markers.
plot_ids = [1, core_ids];
fig7 = figure('Color','w');
tl7 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% left panel: zoom range
ax1 = nexttile; hold(ax1, 'on');
for ii = 1:numel(plot_ids)
    ic = plot_ids(ii);
    D = cases(ic).res.output_data.diam_i(:);
    N = cases(ic).res.output_data.nspec_i(end, :)';
    loglog(ax1, D, N, '-', 'LineWidth', 1.5, 'Color', cols(ic,:), ...
        'DisplayName', case_label(cases(ic)));
end
for ii = 1:numel(core_ids)
    ic = core_ids(ii);
    dmark = cases(ic).dmax_cm_img;
    if isfinite(dmark) && dmark >= plot_range_cm(1) && dmark <= plot_range_cm(2)
        xline(ax1, dmark, '--', 'Color', cols(ic,:), 'LineWidth', 1.0, ...
            'HandleVisibility','off');
    end
end
grid(ax1, 'on');
xlabel(ax1, 'Image diameter (cm)');
ylabel(ax1, 'Number spectrum N');
title(ax1, 'Zoom view');
set(ax1, 'XScale', 'log', 'YScale', 'log');
xlim(ax1, plot_range_cm);

% right panel: wide range to show mapped Dmax markers
dmarks = arrayfun(@(c) c.dmax_cm_img, cases(core_ids));
dmarks = dmarks(isfinite(dmarks) & dmarks > 0);
Dmax_grid = max(cell2mat(arrayfun(@(ic) cases(ic).res.output_data.diam_i(:), plot_ids, 'UniformOutput', false)));
if isempty(dmarks)
    xmax_wide = max(plot_range_cm(2), Dmax_grid);
else
    xmax_wide = max([plot_range_cm(2), Dmax_grid, 1.2 * max(dmarks)]);
end
xrange_wide = [plot_range_cm(1), xmax_wide];

ax2 = nexttile; hold(ax2, 'on');
for ii = 1:numel(plot_ids)
    ic = plot_ids(ii);
    D = cases(ic).res.output_data.diam_i(:);
    N = cases(ic).res.output_data.nspec_i(end, :)';
    loglog(ax2, D, N, '-', 'LineWidth', 1.5, 'Color', cols(ic,:), ...
        'DisplayName', case_label(cases(ic)));
end
for ii = 1:numel(core_ids)
    ic = core_ids(ii);
    dmark = cases(ic).dmax_cm_img;
    if isfinite(dmark) && dmark >= xrange_wide(1) && dmark <= xrange_wide(2)
        xline(ax2, dmark, '--', 'Color', cols(ic,:), 'LineWidth', 1.0, ...
            'HandleVisibility','off');
    end
end
grid(ax2, 'on');
xlabel(ax2, 'Image diameter (cm)');
ylabel(ax2, 'Number spectrum N');
title(ax2, 'Wide view (mapped markers visible)');
set(ax2, 'XScale', 'log', 'YScale', 'log');
xlim(ax2, xrange_wide);
legend(ax2, 'Location', 'best');
title(tl7, 'Final spectra with D_{max} markers (mapped to image diameter, dashed)');
saveas(fig7, fullfile(fig_dir, 'powerlaw_frag_levels_dmax_markers.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);

% ---------- text summary ----------
fprintf('\n=== Final summary ===\n');
for ic = 1:numel(cases)
    fit = cases(ic).fit;
    t = cases(ic).res.time(:);
    mask = t >= steady_t0;
    b_med = median(cases(ic).b_t(mask), 'omitnan');
    b_rng = [min(cases(ic).b_t(mask)), max(cases(ic).b_t(mask))];
    curv = max(abs(fit.resid));
    if isnan(cases(ic).eps)
        fprintf('%s | R2=%.4f | b(final)=%.4f | b(steady med)=%.4f | b-range=[%.4f %.4f] | max|resid|=%.4f | tail(last/prev)=%.3f\n', ...
            case_label(cases(ic)), fit.r2, fit.b, b_med, b_rng(1), b_rng(2), curv, ...
            cases(ic).tail_ratio_last);
    else
        fprintf('%s | eps=%0.0e | Dmax(vol)=%.3f cm | Dmax(img)=%.3f cm | R2=%.4f | b(final)=%.4f | b(steady med)=%.4f | b-range=[%.4f %.4f] | max|resid|=%.4f | tail(last/prev)=%.3f\n', ...
            case_label(cases(ic)), cases(ic).eps, cases(ic).dmax_cm_vol, cases(ic).dmax_cm_img, ...
            fit.r2, fit.b, b_med, b_rng(1), b_rng(2), curv, ...
            cases(ic).tail_ratio_last);
    end
end

if ~isempty(core_ids)
    fprintf('\n=== Core eps slope table ===\n');
    fprintf('eps\tDmax_vol(cm)\tDmax_img(cm)\tb(final)\tb(steady median)\n');
    for ii = 1:numel(core_ids)
        ic = core_ids(ii);
        tt = cases(ic).res.time(:);
        mask = tt >= steady_t0;
        b_med = median(cases(ic).b_t(mask), 'omitnan');
        fprintf('%0.0e\t%.3f\t%.3f\t%.4f\t%.4f\n', ...
            cases(ic).eps, cases(ic).dmax_cm_vol, cases(ic).dmax_cm_img, ...
            cases(ic).fit.b, b_med);
    end
end

% ---------- local helpers ----------
function txt = case_label(c)
if isnan(c.eps)
    txt = 'no frag';
else
    if isfield(c, 'is_core') && ~c.is_core
        txt = sprintf('frag eps=%0.0e (extra)', c.eps);
    else
        txt = sprintf('frag eps=%0.0e', c.eps);
    end
end
end

function fit = fit_powerlaw_loglog(D, N, fit_rng)
fit = struct('ok',false,'a',nan,'b',nan,'slope',nan,'r2',nan, ...
    'Dsel',[],'resid',[],'Dfit',[]);

ok = isfinite(D) & isfinite(N) & (D > 0) & (N > 0) & ...
    (D >= fit_rng(1)) & (D <= fit_rng(2));
if nnz(ok) < 4
    return
end

Dsel = D(ok);
Nsel = N(ok);
x = log10(Dsel);
y = log10(Nsel);

p = polyfit(x, y, 1);
yhat = polyval(p, x);

ss_res = sum((y - yhat).^2);
ss_tot = sum((y - mean(y)).^2);
if ss_tot <= 0
    r2 = NaN;
else
    r2 = 1 - ss_res/ss_tot;
end

a = 10^(p(2));
b = -p(1);

fit.ok = true;
fit.a = a;
fit.b = b;
fit.slope = p(1);
fit.r2 = r2;
fit.Dsel = Dsel;
fit.resid = y - yhat;
fit.Dfit = logspace(log10(min(Dsel)), log10(max(Dsel)), 60);
end

function [edges, dD] = bin_edges_from_centers(D)
% Build bin edges from center diameters (log-midpoint style).
D = D(:);
n = numel(D);
edges = zeros(n+1,1);
if n < 2
    edges(1) = D(1) / 1.2;
    edges(2) = D(1) * 1.2;
    dD = diff(edges);
    return
end

edges(2:n) = sqrt(D(1:n-1) .* D(2:n));
edges(1) = D(1)^2 / edges(2);
edges(n+1) = D(n)^2 / edges(n);
dD = diff(edges);
dD(dD <= 0) = eps;
end

function k = poisson_sample(lambda)
% Simple Poisson sampling with fallback.
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
        % normal approx for faster sampling at high lambda
        k(i) = max(0, round(L + sqrt(L) * randn()));
    end
end
end
