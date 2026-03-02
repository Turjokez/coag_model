% run_uvp_sigma_volume_test.m
% Final UVP-style check:
% - only no-frag, eps=1e-8,1e-6,1e-4
% - UVP sample volume test: 0.1, 1, 10 L
% - sigma residual stats and bins above 2 sigma

clear; close all; clc;
setup_paths
rng(1) % fixed random seed for repeatable Poisson sample

% ---------- output folder ----------
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
fig_dir = fullfile(proj_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% ---------- settings ----------
fit_range_cm = [0.01, 0.6];
plot_range_cm = [0.005, 2.0];
eps_core = [1e-8, 1e-6, 1e-4];
vol_list_L = [0.1, 1, 10];
sinking_law = 'kriest_8';
best_scale = 0.8; % from previous run

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
base.sinking_scale = best_scale;

% ---------- run 4 cases ----------
cases = struct([]);

% no frag
cfg = base.copy();
cfg.enable_disagg = false;
cfg.disagg_mode = 'legacy';
sim = CoagulationSimulation(cfg);
res = sim.run();
cases(1).name = 'no_frag';
cases(1).eps = NaN;
cases(1).res = res;

% frag cases
for j = 1:numel(eps_core)
    cfg = base.copy();
    cfg.enable_disagg = true;
    cfg.disagg_mode = 'operator_split';
    cfg.disagg_epsilon = eps_core(j);
    cfg.disagg_outer_dt = 1/24;
    cfg.disagg_dmax_cm = [];

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    cases(j+1).name = sprintf('frag_eps_%0.0e', eps_core(j));
    cases(j+1).eps = eps_core(j);
    cases(j+1).res = res;
end

% ---------- Figure 1: final model spectra ----------
fig1 = figure('Color','w'); hold on;
cols = lines(numel(cases));
for ic = 1:numel(cases)
    D = cases(ic).res.output_data.diam_i(:);
    N = cases(ic).res.output_data.nspec_i(end, :)' * 1e3; % #/L/cm
    loglog(D, N, 'LineWidth', 1.5, 'Color', cols(ic,:), ...
        'DisplayName', case_label(cases(ic)));
end
grid on;
xlabel('Image diameter (cm)');
ylabel('N (#/L/cm)');
title('Final model spectra (no frag + 3 frag levels)');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(plot_range_cm);
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'uvp_sigma_final_model_spectra.png'));

% ---------- sigma table ----------
sigma_report = struct([]);
krep = 20; % repeats for robust counts

fprintf('\n=== Sigma residual summary (|z| > 2 bins in fit range) ===\n');
fprintf('Each value is median across %d Poisson repeats.\n', krep);
fprintf('Case\\Volume(L)\t0.1\t1\t10\n');

for ic = 1:numel(cases)
    line_vals = nan(1, numel(vol_list_L));
    for iv = 1:numel(vol_list_L)
        V = vol_list_L(iv);
        n_over2 = nan(krep,1);

        for r = 1:krep
            stats = sigma_stats_from_case(cases(ic), V, fit_range_cm);
            n_over2(r) = stats.n_over2;

            % Keep one representative set for plot at V=1 L only
            if (V == 1) && (r == 1)
                sigma_report(ic).V1 = stats;
            end
        end

        sigma_report(ic).name = cases(ic).name;
        sigma_report(ic).V(iv).V = V;
        sigma_report(ic).V(iv).n_over2_med = median(n_over2, 'omitnan');
        sigma_report(ic).V(iv).n_over2_mean = mean(n_over2, 'omitnan');
        line_vals(iv) = sigma_report(ic).V(iv).n_over2_med;
    end

    fprintf('%s\t%.1f\t%.1f\t%.1f\n', case_label(cases(ic)), ...
        line_vals(1), line_vals(2), line_vals(3));
end

% ---------- Figure 2: sampled points + one-particle line at V=1 L ----------
fig2 = figure('Color','w');
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
for ic = 1:numel(cases)
    s = sigma_report(ic).V1;
    D = s.D;
    loglog(D, s.N_model_L, '-', 'Color', cols(ic,:), 'LineWidth', 1.2, ...
        'DisplayName', [case_label(cases(ic)) ' model']);
    loglog(D(s.ok), s.N_obs_L(s.ok), 'o', 'Color', cols(ic,:), 'MarkerSize', 4, ...
        'DisplayName', [case_label(cases(ic)) ' sampled']);
    loglog(D, s.one_particle, '--', 'Color', cols(ic,:), 'LineWidth', 1.0, ...
        'HandleVisibility','off');
end
grid on;
xlabel('Image diameter (cm)');
ylabel('N (#/L/cm)');
title('UVP-like sampled points (V = 1 L)');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(plot_range_cm);
legend('Location','eastoutside');

nexttile; hold on;
for ic = 1:numel(cases)
    s = sigma_report(ic).V1;
    if isempty(s.logD_fit)
        continue
    end
    plot(s.logD_fit, s.z, '-o', 'Color', cols(ic,:), 'LineWidth', 1.1, 'MarkerSize', 4, ...
        'DisplayName', case_label(cases(ic)));
end
yline(2, 'k--', 'LineWidth', 1.0, 'DisplayName','+2 sigma');
yline(-2, 'k--', 'LineWidth', 1.0, 'HandleVisibility','off');
yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
grid on;
xlabel('log10(D)');
ylabel('z = residual / sigma');
title('Sigma residuals (fit window, V = 1 L)');
legend('Location','eastoutside');

title(tl, 'UVP sigma test (no frag + 3 frag levels)');
saveas(fig2, fullfile(fig_dir, 'uvp_sigma_sampled_and_z_V1L.png'));

% ---------- Figure 3: n(|z|>2) vs sample volume ----------
fig3 = figure('Color','w'); hold on;
for ic = 1:numel(cases)
    y = nan(1,numel(vol_list_L));
    for iv = 1:numel(vol_list_L)
        y(iv) = sigma_report(ic).V(iv).n_over2_med;
    end
    semilogx(vol_list_L, y, '-o', 'LineWidth', 1.4, 'Color', cols(ic,:), ...
        'DisplayName', case_label(cases(ic)));
end
grid on;
xlabel('Sample volume (L)');
ylabel('Median # bins with |z| > 2');
title('Detectability vs UVP sample volume');
legend('Location','best');
saveas(fig3, fullfile(fig_dir, 'uvp_sigma_over2_vs_volume.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);

% ---------- helper: sigma stats ----------
function out = sigma_stats_from_case(c, V_L, fit_rng)
    D = c.res.output_data.diam_i(:);
    N_model_L = c.res.output_data.nspec_i(end, :)' * 1e3; % #/L/cm
    [~, dD] = bin_edges_from_centers(D);

    lambda = N_model_L .* dD .* V_L;   % expected counts in bin
    counts = poisson_sample(lambda);
    N_obs_L = counts ./ max(V_L * dD, eps);
    one_particle = 1 ./ max(V_L * dD, eps);

    % fit on sampled points where counts>0 in fit range
    ok = counts > 0 & D >= fit_rng(1) & D <= fit_rng(2);
    fit = fit_powerlaw_loglog(D(ok), N_obs_L(ok), fit_rng);

    z = [];
    logD_fit = [];
    n_over2 = NaN;
    if fit.ok
        x = log10(fit.Dsel);
        y = log10(fit.a * fit.Dsel.^(-fit.b));
        yobs = log10(fit.Nsel);

        % sigma in log-space from Poisson counts: sigma_log10 ~= 1/(ln(10)*sqrt(k))
        kfit = counts(ok);
        kfit = max(kfit, 1);
        sigma_log = 1 ./ (log(10) * sqrt(kfit));
        z = (yobs - y) ./ sigma_log;
        logD_fit = x;
        n_over2 = sum(abs(z) > 2);
    end

    out = struct();
    out.D = D;
    out.dD = dD;
    out.N_model_L = N_model_L;
    out.N_obs_L = N_obs_L;
    out.one_particle = one_particle;
    out.counts = counts;
    out.ok = ok;
    out.fit = fit;
    out.z = z;
    out.logD_fit = logD_fit;
    out.n_over2 = n_over2;
end

% ---------- local helpers ----------
function txt = case_label(c)
if isnan(c.eps)
    txt = 'no frag';
else
    txt = sprintf('frag eps=%0.0e', c.eps);
end
end

function fit = fit_powerlaw_loglog(D, N, fit_rng)
fit = struct('ok',false,'a',nan,'b',nan,'slope',nan,'r2',nan, ...
    'Dsel',[],'Nsel',[],'resid',[],'Dfit',[]);
if isempty(D) || isempty(N)
    return
end
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

fit.ok = true;
fit.a = 10^(p(2));
fit.b = -p(1);
fit.slope = p(1);
fit.r2 = r2;
fit.Dsel = Dsel;
fit.Nsel = Nsel;
fit.resid = y - yhat;
fit.Dfit = logspace(log10(min(Dsel)), log10(max(Dsel)), 60);
end

function [edges, dD] = bin_edges_from_centers(D)
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
