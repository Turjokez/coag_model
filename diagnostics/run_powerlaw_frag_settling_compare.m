% run_powerlaw_frag_settling_compare.m
% Compare power-law + fragmentation results for two settling laws.
% Also track average sinking speed over time in each run.

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
fit_range_cm = [0.01, 0.6];
plot_range_cm = [0.005, 2.0];
eps_core = [1e-8, 1e-6, 1e-4];
scale_scan = [0.8, 1.0, 1.2];
steady_t0 = 20;
sinking_laws = {'kriest_8', 'current'}; % current = older law

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
base.sinking_size = 'image';

all_runs = struct([]);

for ilaw = 1:numel(sinking_laws)
    law = sinking_laws{ilaw};
    fprintf('\n==============================\n');
    fprintf('Settling law: %s\n', law);
    fprintf('==============================\n');

    base_law = base.copy();
    base_law.sinking_law = law;

    % ---- Step 1: pick best no-frag baseline for this settling law ----
    scan = struct([]);
    for i = 1:numel(scale_scan)
        cfg = base_law.copy();
        cfg.sinking_scale = scale_scan(i);
        cfg.enable_disagg = false;
        cfg.disagg_mode = 'legacy';

        sim = CoagulationSimulation(cfg);
        res = sim.run();

        D = res.output_data.diam_i(:);
        N = res.output_data.nspec_i(end, :)';
        fit = fit_powerlaw_loglog(D, N, fit_range_cm);

        scan(i).cfg = cfg;
        scan(i).res = res;
        scan(i).fit = fit;
        scan(i).scale = scale_scan(i);

        fprintf('scale=%.2f | R2=%.4f | b=%.4f\n', scale_scan(i), fit.r2, fit.b);
    end

    [~, best_i] = max(arrayfun(@(s) s.fit.r2, scan));
    best_scale = scan(best_i).scale;
    res_no_frag = scan(best_i).res;
    fit_no_frag = scan(best_i).fit;
    fprintf('Selected baseline for %s: scale=%.2f (R2=%.4f)\n', ...
        law, best_scale, fit_no_frag.r2);

    % ---- Step 2/3: no-frag + three frag levels ----
    cases = struct([]);

    % no frag
    cases(1).name = 'no_frag';
    cases(1).eps = NaN;
    cases(1).cfg = scan(best_i).cfg;
    cases(1).res = res_no_frag;
    cases(1).fit = fit_no_frag;

    % frag levels
    for j = 1:numel(eps_core)
        idx = j + 1;
        cfg = base_law.copy();
        cfg.sinking_scale = best_scale;
        cfg.enable_disagg = true;
        cfg.disagg_mode = 'operator_split';
        cfg.disagg_epsilon = eps_core(j);
        cfg.disagg_outer_dt = 1/24;
        cfg.disagg_dmax_cm = [];

        sim = CoagulationSimulation(cfg);
        res = sim.run();

        D = res.output_data.diam_i(:);
        N = res.output_data.nspec_i(end, :)';
        fit = fit_powerlaw_loglog(D, N, fit_range_cm);

        cases(idx).name = sprintf('frag_eps_%0.0e', eps_core(j));
        cases(idx).eps = eps_core(j);
        cases(idx).cfg = cfg;
        cases(idx).res = res;
        cases(idx).fit = fit;

        fprintf('eps=%0.0e | R2=%.4f | b=%.4f\n', eps_core(j), fit.r2, fit.b);
    end

    % ---- diagnostics over time (b + average sinking speed) ----
    for ic = 1:numel(cases)
        t = cases(ic).res.time(:);
        D = cases(ic).res.output_data.diam_i(:);
        Ns = cases(ic).res.output_data.nspec_i;      % number spectrum (time x bins)
        Y = max(cases(ic).res.concentrations, 0);    % state (time x bins)
        w = cases(ic).res.output_data.set_vel(:)';   % m/day (1 x bins)

        b_t = nan(numel(t),1);
        for k = 1:numel(t)
            fit_k = fit_powerlaw_loglog(D, Ns(k,:)', fit_range_cm);
            b_t(k) = fit_k.b;
        end
        cases(ic).b_t = b_t;

        % Number-weighted mean settling speed
        den_num = sum(max(Ns, 0), 2);
        num_num = max(Ns, 0) * w';
        wbar_num = num_num ./ max(den_num, eps);

        % Mass/state-weighted mean settling speed
        den_mass = sum(Y, 2);
        num_mass = Y * w';
        wbar_mass = num_mass ./ max(den_mass, eps);

        cases(ic).wbar_num_t = wbar_num;
        cases(ic).wbar_mass_t = wbar_mass;

        mask = t >= steady_t0;
        cases(ic).b_med_steady = median(b_t(mask), 'omitnan');
        cases(ic).wbar_num_steady = median(wbar_num(mask), 'omitnan');
        cases(ic).wbar_mass_steady = median(wbar_mass(mask), 'omitnan');
    end

    all_runs(ilaw).law = law;
    all_runs(ilaw).best_scale = best_scale;
    all_runs(ilaw).cases = cases;
end

% ---------- Figure 1: final spectra (two-law comparison) ----------
fig1 = figure('Color','w');
tl1 = tiledlayout(1, numel(all_runs), 'TileSpacing', 'compact', 'Padding', 'compact');
case_cols = [0.10 0.45 0.80; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.86];
for ilaw = 1:numel(all_runs)
    nexttile; hold on;
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        D = cases(ic).res.output_data.diam_i(:);
        N = cases(ic).res.output_data.nspec_i(end, :)';
        loglog(D, N, '-', 'Color', case_cols(ic,:), 'LineWidth', 1.5, ...
            'DisplayName', case_label(cases(ic)));
    end
    grid on;
    xlabel('Image diameter (cm)');
    ylabel('Number spectrum N');
    title(sprintf('%s (scale=%.2f)', all_runs(ilaw).law, all_runs(ilaw).best_scale), ...
        'Interpreter','none');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim(plot_range_cm);
    legend('Location', 'best');
end
title(tl1, 'Final spectra: no-frag + 3 frag levels (two settling laws)');
saveas(fig1, fullfile(fig_dir, 'settling_compare_final_spectra.png'));

% ---------- Figure 2: slope vs epsilon (two laws) ----------
fig2 = figure('Color','w');
tl2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% b(final)
nexttile; hold on;
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    epsv = [cases(2:end).eps];
    bfin = arrayfun(@(c) c.fit.b, cases(2:end));
    plot(epsv, bfin, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('%s frag', all_runs(ilaw).law));
    yline(cases(1).fit.b, '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s no-frag b', all_runs(ilaw).law));
end
set(gca, 'XScale', 'log');
grid on;
xlabel('Fragmentation epsilon');
ylabel('b (final)');
title('Final slope vs epsilon');
legend('Location','best');

% b(steady median)
nexttile; hold on;
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    epsv = [cases(2:end).eps];
    bmed = arrayfun(@(c) c.b_med_steady, cases(2:end));
    plot(epsv, bmed, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('%s frag', all_runs(ilaw).law));
    yline(cases(1).b_med_steady, '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s no-frag med b', all_runs(ilaw).law));
end
set(gca, 'XScale', 'log');
grid on;
xlabel('Fragmentation epsilon');
ylabel(sprintf('b median (t >= %g d)', steady_t0));
title('Steady slope vs epsilon');
legend('Location','best');

title(tl2, 'Slope comparison between settling laws');
saveas(fig2, fullfile(fig_dir, 'settling_compare_slope_vs_epsilon.png'));

% ---------- Figure 3: average sinking speed over time ----------
fig3 = figure('Color','w');
tl3 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
law_styles = {'-', '--'};

% number-weighted
nexttile; hold on;
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        t = cases(ic).res.time(:);
        plot(t, cases(ic).wbar_num_t, 'LineStyle', law_styles{ilaw}, ...
            'Color', case_cols(ic,:), 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s | %s', all_runs(ilaw).law, case_label(cases(ic))));
    end
end
grid on;
xlabel('Time (days)');
ylabel('Avg sinking speed (m/day)');
title('Number-weighted mean settling speed');
legend('Location','eastoutside');

% mass-weighted
nexttile; hold on;
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        t = cases(ic).res.time(:);
        plot(t, cases(ic).wbar_mass_t, 'LineStyle', law_styles{ilaw}, ...
            'Color', case_cols(ic,:), 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s | %s', all_runs(ilaw).law, case_label(cases(ic))));
    end
end
grid on;
xlabel('Time (days)');
ylabel('Avg sinking speed (m/day)');
title('Mass-weighted mean settling speed');
legend('Location','eastoutside');

title(tl3, 'Average settling speed over time (all runs)');
saveas(fig3, fullfile(fig_dir, 'settling_compare_avg_sinking_speed_timeseries.png'));

% ---------- Figure 4: steady mean sinking speed summary ----------
fig4 = figure('Color','w');
tl4 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
case_names = {'no frag', 'eps=1e-8', 'eps=1e-6', 'eps=1e-4'};
x = 1:numel(case_names);

% number-weighted steady summary
nexttile; hold on;
Ynum = nan(numel(all_runs), numel(case_names));
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        Ynum(ilaw, ic) = cases(ic).wbar_num_steady;
    end
end
b1 = bar(x - 0.18, Ynum(1,:), 0.35, 'FaceColor', [0.10 0.45 0.80]);
b2 = bar(x + 0.18, Ynum(2,:), 0.35, 'FaceColor', [0.85 0.33 0.10]);
grid on;
set(gca, 'XTick', x, 'XTickLabel', case_names);
ylabel(sprintf('Median avg speed (m/day), t >= %g d', steady_t0));
title('Number-weighted steady mean speed');
legend([b1 b2], {all_runs(1).law, all_runs(2).law}, 'Location', 'best');

% mass-weighted steady summary
nexttile; hold on;
Ymass = nan(numel(all_runs), numel(case_names));
for ilaw = 1:numel(all_runs)
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        Ymass(ilaw, ic) = cases(ic).wbar_mass_steady;
    end
end
b3 = bar(x - 0.18, Ymass(1,:), 0.35, 'FaceColor', [0.10 0.45 0.80]);
b4 = bar(x + 0.18, Ymass(2,:), 0.35, 'FaceColor', [0.85 0.33 0.10]);
grid on;
set(gca, 'XTick', x, 'XTickLabel', case_names);
ylabel(sprintf('Median avg speed (m/day), t >= %g d', steady_t0));
title('Mass-weighted steady mean speed');
legend([b3 b4], {all_runs(1).law, all_runs(2).law}, 'Location', 'best');

title(tl4, 'Steady average sinking speed summary');
saveas(fig4, fullfile(fig_dir, 'settling_compare_avg_sinking_speed_steady.png'));

fprintf('\nSaved figures to: %s\n', fig_dir);

% ---------- text summary ----------
fprintf('\n=== Settling-law comparison summary ===\n');
for ilaw = 1:numel(all_runs)
    fprintf('\nLaw: %s | best scale=%.2f\n', all_runs(ilaw).law, all_runs(ilaw).best_scale);
    cases = all_runs(ilaw).cases;
    for ic = 1:numel(cases)
        if isnan(cases(ic).eps)
            fprintf('  %-14s | b(final)=%.4f | b(steady)=%.4f | w_num(steady)=%.3f | w_mass(steady)=%.3f m/day\n', ...
                case_label(cases(ic)), cases(ic).fit.b, cases(ic).b_med_steady, ...
                cases(ic).wbar_num_steady, cases(ic).wbar_mass_steady);
        else
            fprintf('  %-14s | eps=%0.0e | b(final)=%.4f | b(steady)=%.4f | w_num(steady)=%.3f | w_mass(steady)=%.3f m/day\n', ...
                case_label(cases(ic)), cases(ic).eps, cases(ic).fit.b, cases(ic).b_med_steady, ...
                cases(ic).wbar_num_steady, cases(ic).wbar_mass_steady);
        end
    end
end

% ---------- local helpers ----------
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

fit.ok = true;
fit.a = 10^(p(2));
fit.b = -p(1);
fit.slope = p(1);
fit.r2 = r2;
fit.Dsel = Dsel;
fit.resid = y - yhat;
fit.Dfit = logspace(log10(min(Dsel)), log10(max(Dsel)), 60);
end

function txt = case_label(c)
if isnan(c.eps)
    txt = 'no frag';
else
    txt = sprintf('frag %0.0e', c.eps);
end
end
