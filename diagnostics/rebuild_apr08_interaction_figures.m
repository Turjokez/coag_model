function rebuild_apr08_interaction_figures
% rebuild_apr08_interaction_figures
% Short note:
% 1. rebuild the Apr 08 interaction figures
% 2. keep labels short and easy
% 3. save all figures to docs/figures

clear;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
fig_dir = fullfile(proj_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

d_um = logspace(0, 4, 400);
d_cm = d_um * 1e-4;

plot_exponent_sweep(d_um, d_cm, fig_dir);
plot_repo_settling_laws(d_um, d_cm, fig_dir);
plot_fixed_d1_curves(fig_dir);
plot_named_maps(fig_dir);
plot_ratio_maps(fig_dir);

disp('Saved Apr 08 interaction figures to:');
disp(fig_dir);
end

function plot_exponent_sweep(d_um, d_cm, fig_dir)
b_list = 0.6:0.2:2.0;
cols = parula(numel(b_list));

fig = figure('Color', 'w', 'Position', [80 80 760 480]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(b_list)
    b = b_list(i);
    w = d_cm .^ b;
    loglog(ax, d_um, w, 'LineWidth', 2.0, 'Color', cols(i, :), ...
        'DisplayName', sprintf('b = %.1f', b));
end

set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [min(d_um) max(d_um)]);
grid(ax, 'on');
xlabel(ax, 'Particle diameter, d (um)');
ylabel(ax, 'Relative sinking speed, w');
legend(ax, 'Location', 'southeast', 'NumColumns', 2, 'Box', 'off');
plain_style(fig);
save_png(fig, fullfile(fig_dir, 'apr08_sinking_speed_exponent_sweep.png'));
end

function plot_repo_settling_laws(d_um, d_cm, fig_dir)
laws = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
cols = lines(numel(laws));

fig = figure('Color', 'w', 'Position', [80 80 740 460]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(laws)
    w = named_speed_mday(d_cm, laws{i});
    loglog(ax, d_um, w, 'LineWidth', 2.0, 'Color', cols(i, :), ...
        'DisplayName', laws{i});
end

set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [min(d_um) max(d_um)]);
grid(ax, 'on');
xlabel(ax, 'Particle diameter, d (um)');
ylabel(ax, 'Sinking speed (m day^{-1})');
legend(ax, 'Location', 'northwest', 'Box', 'off');
plain_style(fig);
save_png(fig, fullfile(fig_dir, 'apr08_repo_settling_speed_laws.png'));
end

function plot_fixed_d1_curves(fig_dir)
d2_um = logspace(0, 4, 400);
d2_cm = d2_um * 1e-4;
b_list = 0.6:0.2:2.0;
cols = parula(numel(b_list));

cases = {
    1, '1um'
    500, '500um'
    1000, '1mm'
    };

for k = 1:size(cases, 1)
    d1_um = cases{k, 1};
    tag = cases{k, 2};
    d1_cm = d1_um * 1e-4;

    fig = figure('Color', 'w', 'Position', [80 80 760 480]);
    ax = axes(fig);
    hold(ax, 'on');

    for i = 1:numel(b_list)
        b = b_list(i);
        w1 = d1_cm .^ b;
        w2 = d2_cm .^ b;
        beta = beta_diff_sed(d1_cm, d2_cm, w1, w2);
        beta(beta <= 0) = NaN;
        loglog(ax, d2_um, beta, 'LineWidth', 2.0, 'Color', cols(i, :), ...
            'DisplayName', sprintf('b = %.1f', b));
    end

    set(ax, 'XScale', 'log', 'YScale', 'log');
    xlim(ax, [min(d2_um) max(d2_um)]);
    grid(ax, 'on');
    title(ax, sprintf('Fixed particle, d1 = %g um', d1_um), 'FontWeight', 'normal');
    xlabel(ax, 'Partner particle size, d2 (um)');
    ylabel(ax, 'Differential-settling beta (comparison units)');
    legend(ax, 'Location', 'northwest', 'NumColumns', 2, 'Box', 'off');
    plain_style(fig);
    save_png(fig, fullfile(fig_dir, ['apr08_beta_fixed_d1_' tag '.png']));
end
end

function plot_named_maps(fig_dir)
d_um = logspace(0, 4, 220);
d_cm = d_um * 1e-4;
[D1_cm, D2_cm] = meshgrid(d_cm, d_cm);
[D1_um, D2_um] = meshgrid(d_um, d_um);

laws = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
beta_maps = cell(1, numel(laws));
rate_maps = cell(1, numel(laws));

C = powerlaw_concentration(d_cm, 1e3, -2.5);
[C1, C2] = meshgrid(C, C);

all_beta = [];
all_rate = [];

for i = 1:numel(laws)
    w = named_speed_cms(d_cm, laws{i});
    [W1, W2] = meshgrid(w, w);
    beta = beta_diff_sed(D1_cm, D2_cm, W1, W2);
    rate = beta .* C1 .* C2;

    z_beta = log10(max(beta, realmin));
    z_rate = log10(max(rate, realmin));

    beta_maps{i} = z_beta;
    rate_maps{i} = z_rate;

    all_beta = [all_beta; z_beta(isfinite(z_beta))]; %#ok<AGROW>
    all_rate = [all_rate; z_rate(isfinite(z_rate))]; %#ok<AGROW>
end

beta_vmin = prctile(all_beta, 1);
beta_vmax = prctile(all_beta, 99);
rate_vmin = prctile(all_rate, 1);
rate_vmax = prctile(all_rate, 99);

beta_levels = linspace(beta_vmin, beta_vmax, 18);
rate_levels = linspace(rate_vmin, rate_vmax, 18);

fig1 = figure('Color', 'w', 'Position', [80 80 820 720]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_last = [];
for i = 1:numel(laws)
    ax = nexttile(tl1, i);
    ax_last = ax;
    z = min(max(beta_maps{i}, beta_vmin), beta_vmax);
    contourf(ax, D1_um, D2_um, z, beta_levels, 'LineColor', 'none');
    hold(ax, 'on');
    plot(ax, d_um, d_um, 'w-', 'LineWidth', 0.8);
    set(ax, 'XScale', 'log', 'YScale', 'log');
    xlabel(ax, 'd1 (um)');
    ylabel(ax, 'd2 (um)');
    title(ax, laws{i}, 'Interpreter', 'none', 'FontWeight', 'normal');
    grid(ax, 'on');
end
colormap(fig1, parula(256));
cb1 = colorbar(ax_last, 'eastoutside');
cb1.Label.String = 'log10(beta_{DS})';
plain_style(fig1);
save_png(fig1, fullfile(fig_dir, 'apr08_beta_only_repo_laws.png'));

fig2 = figure('Color', 'w', 'Position', [80 80 820 720]);
tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_last = [];
for i = 1:numel(laws)
    ax = nexttile(tl2, i);
    ax_last = ax;
    z = min(max(rate_maps{i}, rate_vmin), rate_vmax);
    contourf(ax, D1_um, D2_um, z, rate_levels, 'LineColor', 'none');
    hold(ax, 'on');
    plot(ax, d_um, d_um, 'Color', [0.45 0.00 0.60], 'LineWidth', 0.8);
    set(ax, 'XScale', 'log', 'YScale', 'log');
    xlabel(ax, 'd1 (um)');
    ylabel(ax, 'd2 (um)');
    title(ax, laws{i}, 'Interpreter', 'none', 'FontWeight', 'normal');
    grid(ax, 'on');
end
colormap(fig2, parula(256));
cb2 = colorbar(ax_last, 'eastoutside');
cb2.Label.String = 'log10(beta_{DS} C1 C2)';
plain_style(fig2);
save_png(fig2, fullfile(fig_dir, 'apr08_beta_c1c2_heatmaps_named_laws.png'));
end

function plot_ratio_maps(fig_dir)
d_um = logspace(0, 4, 220);
d_cm = d_um * 1e-4;
[D1_cm, D2_cm] = meshgrid(d_cm, d_cm);
[D1_um, D2_um] = meshgrid(d_um, d_um);

laws = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
eps_list = [1e-8, 1e-6, 1e-4];

for e = 1:numel(eps_list)
    eps_now = eps_list(e);
    ratio_maps = cell(1, numel(laws));
    all_ratio = [];

    for i = 1:numel(laws)
        w = named_speed_cms(d_cm, laws{i});
        [W1, W2] = meshgrid(w, w);
        beta_ds = beta_diff_sed(D1_cm, D2_cm, W1, W2);
        beta_sh = beta_turb_shear(D1_cm, D2_cm, eps_now);
        ratio = beta_ds ./ max(beta_sh, realmin);
        z = log10(max(ratio, realmin));
        ratio_maps{i} = z;
        all_ratio = [all_ratio; z(isfinite(z))]; %#ok<AGROW>
    end

    q1 = prctile(all_ratio, 1);
    q99 = prctile(all_ratio, 99);
    vmax = max(abs([q1, q99]));
    levels = linspace(-vmax, vmax, 19);

    fig = figure('Color', 'w', 'Position', [80 80 820 720]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_last = [];

    for i = 1:numel(laws)
        ax = nexttile(tl, i);
        ax_last = ax;
        z = min(max(ratio_maps{i}, -vmax), vmax);
        contourf(ax, D1_um, D2_um, z, levels, 'LineColor', 'none');
        hold(ax, 'on');
        contour(ax, D1_um, D2_um, z, [0 0], 'k-', 'LineWidth', 0.8);
        set(ax, 'XScale', 'log', 'YScale', 'log');
        xlabel(ax, 'd1 (um)');
        ylabel(ax, 'd2 (um)');
        title(ax, laws{i}, 'Interpreter', 'none', 'FontWeight', 'normal');
        grid(ax, 'on');
    end

    colormap(fig, div_cmap(256));
    cb = colorbar(ax_last, 'eastoutside');
    cb.Label.String = 'log10(DS / shear)';
    plain_style(fig);
    tag = strrep(sprintf('%.0e', eps_now), 'e-0', 'em0');
    tag = strrep(tag, 'e-', 'em');
    save_png(fig, fullfile(fig_dir, ['apr08_diffsed_vs_shear_ratio_eps_' tag '.png']));
end
end

function w_mday = named_speed_mday(d_cm, law_name)
cfg = local_cfg();
w_mday = named_speed_cms(d_cm, law_name) * cfg.day_to_sec / 100;
end

function w_cms = named_speed_cms(d_cm, law_name)
cfg = local_cfg();
switch lower(string(law_name))
    case "current"
        w_cms = current_speed_cms(d_cm);
    case "kriest_8"
        w_mday = 66 .* d_cm .^ 0.62;
        w_cms = (w_mday * 100) / cfg.day_to_sec;
    case "kriest_9"
        w_mday = 132 .* d_cm .^ 0.62;
        w_cms = (w_mday * 100) / cfg.day_to_sec;
    case "siegel_2025"
        d_mm = d_cm * 10.0;
        w_mday = 20.2 .* d_mm .^ 0.67;
        w_cms = (w_mday * 100) / cfg.day_to_sec;
    otherwise
        error('Unknown law: %s', law_name);
end
end

function w_cms = current_speed_cms(d_cm)
cfg = local_cfg();
r_i = d_cm ./ (2.0 * cfg.r_to_rg);
V = (r_i ./ current_amfrac()) .^ cfg.fr_dim;
r_v = ((0.75 / pi) .* V) .^ (1.0 / 3.0);
w_cms = current_setcon() .* (r_v .^ 3) ./ max(r_i, realmin);
end

function amfrac = current_amfrac()
cfg = local_cfg();
a0 = cfg.d0 / 2.0;
amfrac_temp = (4.0 / 3.0 * pi) ^ (-1.0 / cfg.fr_dim) * a0 ^ (1.0 - 3.0 / cfg.fr_dim);
amfrac = amfrac_temp * sqrt(0.6);
end

function setcon = current_setcon()
cfg = local_cfg();
del_rho = (4.5 * 2.48) * cfg.kvisc * cfg.rho_fl / cfg.g * (cfg.d0 / 2.0) ^ (-0.83);
setcon = (2.0 / 9.0) * del_rho / cfg.rho_fl * cfg.g / cfg.kvisc;
end

function beta = beta_diff_sed(d1_cm, d2_cm, w1_cms, w2_cms)
geom = (pi / 4.0) .* (d1_cm + d2_cm) .^ 2;
beta = geom .* abs(w1_cms - w2_cms);
end

function beta = beta_turb_shear(d1_cm, d2_cm, eps_mks)
cfg = local_cfg();
r1 = d1_cm ./ (2.0 * cfg.r_to_rg);
r2 = d2_cm ./ (2.0 * cfg.r_to_rg);
p = min(r1, r2) ./ max(r1, r2);
p1 = 1.0 + p;
eff = 1.0 - (1.0 + 5.0 .* p + 2.5 .* p .* p) ./ (p1 .^ 5);
rg = (r1 + r2) .* cfg.r_to_rg;
shape = sqrt(8.0 * pi / 15.0) .* eff .* rg .^ 3;
eps_cgs = eps_mks * 1e4;
gamma = sqrt(eps_cgs / cfg.kvisc);
beta = shape .* gamma;
end

function C = powerlaw_concentration(d_cm, amp, expo)
C = amp .* d_cm .^ expo;
end

function plain_style(fig)
ax = findall(fig, 'Type', 'axes');
for k = 1:numel(ax)
    set(ax(k), 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 0.8);
    set(ax(k).Title, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');
    set(ax(k).XLabel, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');
    set(ax(k).YLabel, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');
end

lg = findall(fig, 'Type', 'legend');
for k = 1:numel(lg)
    set(lg(k), 'FontName', 'Arial', 'FontSize', 9);
end

cb = findall(fig, 'Type', 'ColorBar');
for k = 1:numel(cb)
    set(cb(k), 'FontName', 'Arial', 'FontSize', 9);
    set(cb(k).Label, 'FontName', 'Arial', 'FontSize', 9);
end
end

function cmap = div_cmap(n)
if nargin < 1
    n = 256;
end
n1 = floor(n / 2);
n2 = n - n1;
neg = [linspace(0.25, 0.85, n1)', linspace(0.35, 0.92, n1)', ones(n1, 1)];
pos = [ones(n2, 1), linspace(0.92, 0.20, n2)', linspace(0.92, 0.20, n2)'];
cmap = [neg; pos];
end

function save_png(fig, out_path)
set(fig, 'PaperPositionMode', 'auto');
print(fig, out_path, '-dpng', '-r220');
close(fig);
end

function cfg = local_cfg()
cfg = struct();
cfg.r_to_rg = 1.6;
cfg.rho_fl = 1.0275;
cfg.kvisc = 0.01;
cfg.g = 980.0;
cfg.day_to_sec = 86400.0;
cfg.d0 = 20e-4;
cfg.fr_dim = 2.33;
end
