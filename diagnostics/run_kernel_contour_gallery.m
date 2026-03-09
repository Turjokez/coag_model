% run_kernel_contour_gallery.m
% Kernel figures from Adrian handnote.

clear; close all; clc;
script_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_dir);
addpath(proj_root);
run(fullfile(proj_root, 'setup_paths.m'));
set(groot, 'defaultAxesFontName', 'Calibri');
set(groot, 'defaultTextFontName', 'Calibri');
set(groot, 'defaultLegendFontName', 'Calibri');
set(groot, 'defaultColorbarFontName', 'Calibri');
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
set(groot, 'defaultLegendFontWeight', 'normal');
set(groot, 'defaultColorbarFontWeight', 'normal');

fig_dir = fullfile(proj_root, 'output', 'figures');
docs_fig_dir = fullfile(proj_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(docs_fig_dir, 'dir'), mkdir(docs_fig_dir); end

cfg = SimulationConfig();
cfg.sinking_size = 'image';
grid = cfg.derive();

size_um = logspace(0, log10(5000), 180);
fixed_um = [1, 50, 100, 500, 1000, 5000];
temp_C = linspace(0, 30, 140);
eps_mks = logspace(-10, -4, 140); % W/kg = m^2/s^3

% Main handnote compare
law_list = {'current', 'kriest_8', 'laurenceau_cornec', 'white'};
law_labels = {'orig/current V_s', 'kriest V_s', 'Laurenceau-Cornec case', 'White drag case'};

% Extra repo compare
repo_law_list = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
repo_law_labels = {'orig/current V_s', 'kriest_8 V_s', 'kriest_9 V_s', 'siegel_2025 V_s'};

% -------- Brownian: temp vs partner size for fixed particle size --------
brown_z = cell(numel(fixed_um), 1);
brown_all = [];
for k = 1:numel(fixed_um)
    [D2, TT] = meshgrid(size_um, temp_C);
    D1 = fixed_um(k) * ones(size(D2));
    beta_day = brownian_kernel_day(D1, D2, TT, cfg);
    brown_z{k} = log10(max(beta_day, realmin));
    brown_all = [brown_all; brown_z{k}(:)]; %#ok<AGROW>
end
brown_levels = linspace(min(brown_all), max(brown_all), 22);

fig1 = figure('Color', 'w', 'Position', [80 80 1600 850]);
for k = 1:numel(fixed_um)
    [D2, TT] = meshgrid(size_um, temp_C);

    subplot(2, 3, k);
    contourf(D2, TT, brown_z{k}, brown_levels, 'LineColor', 'none');
    colormap(gca, mono_cmap(256, [0.05 0.30 0.85]));
    clim([brown_levels(1), brown_levels(end)]);
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Partner size (um)');
    ylabel('Temperature (^\circC)');
    title(sprintf('Brownian | particle 1 fixed = %g um', fixed_um(k)));
    cb = colorbar;
    cb.Label.String = 'log_{10}(\beta_B) [cm^3 day^{-1}]';
end
h = sgtitle('Brownian coagulation kernel');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig1);
saveas(fig1, fullfile(fig_dir, 'kernel_brownian_temp_size_contours.png'));

% -------- Shear: epsilon vs partner size for fixed particle size --------
shear_z = cell(numel(fixed_um), 1);
shear_all = [];
for k = 1:numel(fixed_um)
    [D2, EE] = meshgrid(size_um, eps_mks);
    D1 = fixed_um(k) * ones(size(D2));
    beta_day = shear_kernel_day(D1, D2, EE, cfg, grid);
    shear_z{k} = log10(max(beta_day, realmin));
    shear_all = [shear_all; shear_z{k}(:)]; %#ok<AGROW>
end
shear_cmax = max(abs(shear_all));
shear_levels = linspace(-shear_cmax, shear_cmax, 23);

fig2 = figure('Color', 'w', 'Position', [80 80 1600 850]);
for k = 1:numel(fixed_um)
    [D2, EE] = meshgrid(size_um, eps_mks);

    subplot(2, 3, k);
    contourf(D2, EE, shear_z{k}, shear_levels, 'LineColor', 'none');
    colormap(gca, div_cmap(256));
    clim([-shear_cmax, shear_cmax]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('Partner size (um)');
    ylabel('\epsilon (W kg^{-1})');
    title(sprintf('Shear | particle 1 fixed = %g um', fixed_um(k)));
    cb = colorbar;
    cb.Label.String = 'log_{10}(\beta_S) [cm^3 day^{-1}]';
end
h = sgtitle('Turbulent shear coagulation kernel');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig2);
saveas(fig2, fullfile(fig_dir, 'kernel_shear_eps_size_contours.png'));

% -------- Differential sedimentation: size-size for handnote laws --------
[D1m, D2m] = meshgrid(size_um, size_um);
ds_main_z = cell(numel(law_list), 1);
ds_main_all = [];
for k = 1:numel(law_list)
    beta_day = ds_kernel_day(D1m, D2m, law_list{k}, cfg, grid);
    ds_main_z{k} = ds_plot_field(beta_day, D1m, D2m, 0.03);
    ds_main_all = [ds_main_all; ds_main_z{k}(:)]; %#ok<AGROW>
end
ds_main_cmax = max(abs(ds_main_all));
ds_main_levels = linspace(-ds_main_cmax, ds_main_cmax, 23);

fig3 = figure('Color', 'w', 'Position', [80 80 1500 1200]);
for k = 1:numel(law_list)
    subplot(2, 2, k);
    contourf(D1m, D2m, ds_main_z{k}, ds_main_levels, 'LineColor', 'none');
    colormap(gca, div_cmap(256));
    clim([-ds_main_cmax, ds_main_cmax]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('Particle 1 size (um)');
    ylabel('Particle 2 size (um)');
    title(sprintf('Differential sedimentation | %s', law_labels{k}));
    cb = colorbar;
    cb.Label.String = 'log_{10}(\beta_{DS}) [cm^3 day^{-1}]';
end
h = sgtitle('Differential sedimentation kernel for handnote settling laws');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig3);
saveas(fig3, fullfile(fig_dir, 'kernel_ds_size_size_settling_laws.png'));

% -------- Extra: settling speed curves --------
fig4 = figure('Color', 'w', 'Position', [100 100 1200 520]);

subplot(1, 2, 1); hold on;
cols = lines(numel(law_list));
for k = 1:numel(law_list)
    w_mday = settling_speed_mday(size_um, law_list{k}, cfg, grid);
    plot(size_um, w_mday, 'LineWidth', 1.6, 'Color', cols(k,:), ...
        'DisplayName', law_labels{k});
end
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Image size (um)');
ylabel('Settling speed (m day^{-1})');
title('Handnote settling speed compare');
legend('Location', 'best');

subplot(1, 2, 2); hold on;
cols = lines(numel(repo_law_list));
for k = 1:numel(repo_law_list)
    w_mday = settling_speed_mday(size_um, repo_law_list{k}, cfg, grid);
    plot(size_um, w_mday, 'LineWidth', 1.6, 'Color', cols(k,:), ...
        'DisplayName', repo_law_labels{k});
end
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Image size (um)');
ylabel('Settling speed (m day^{-1})');
title('Extra repo settling speed compare');
legend('Location', 'best');
h = sgtitle('Settling speed laws');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig4);
saveas(fig4, fullfile(fig_dir, 'kernel_ds_settling_speed_laws.png'));

% -------- Extra: handnote-law ratio to current --------
fig5 = figure('Color', 'w', 'Position', [80 80 1500 420]);
beta_ref = ds_kernel_day(D1m, D2m, 'current', cfg, grid);
extra_laws = law_list(2:end);
extra_labels = law_labels(2:end);
ratio_z = cell(numel(extra_laws), 1);
ratio_all = [];
for k = 1:numel(extra_laws)
    beta_cmp = ds_kernel_day(D1m, D2m, extra_laws{k}, cfg, grid);
    ratio = beta_cmp ./ max(beta_ref, realmin);
    ratio_z{k} = ds_plot_field(ratio, D1m, D2m, 0.03);
    ratio_all = [ratio_all; ratio_z{k}(:)]; %#ok<AGROW>
end
ratio_cmax = max(abs(ratio_all));
ratio_levels = linspace(-ratio_cmax, ratio_cmax, 23);

for k = 1:numel(extra_laws)
    subplot(1, 3, k);
    contourf(D1m, D2m, ratio_z{k}, ratio_levels, 'LineColor', 'none');
    colormap(gca, div_cmap(256));
    clim([-ratio_cmax, ratio_cmax]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('Particle 1 size (um)');
    ylabel('Particle 2 size (um)');
    title(sprintf('%s / current', extra_labels{k}));
    cb = colorbar;
    cb.Label.String = 'log_{10}(ratio)';
end
h = sgtitle('Extra: handnote DS kernel ratio relative to current law');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig5);
saveas(fig5, fullfile(fig_dir, 'kernel_ds_ratio_to_current.png'));

% -------- Extra: DS line cuts for fixed sizes --------
fig6 = figure('Color', 'w', 'Position', [80 80 1500 850]);
cols = lines(numel(law_list));
for k = 1:numel(fixed_um)
    subplot(2, 3, k); hold on;
    y_all = [];
    for j = 1:numel(law_list)
        beta_line = ds_kernel_day(fixed_um(k), size_um, law_list{j}, cfg, grid);
        mask_eq = abs(log10(size_um / fixed_um(k))) < 0.02;
        beta_line(mask_eq) = NaN;
        plot(size_um, beta_line, 'LineWidth', 1.4, 'Color', cols(j,:), ...
            'DisplayName', law_labels{j});
        y_all = [y_all; beta_line(:)]; %#ok<AGROW>
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    y_all = y_all(isfinite(y_all) & y_all > 0);
    if ~isempty(y_all)
        ylim([10^floor(log10(prctile(y_all, 5))), 10^ceil(log10(prctile(y_all, 95)))]);
    end
    grid on;
    xlabel('Partner size (um)');
    ylabel('\beta_{DS} (cm^3 day^{-1})');
    title(sprintf('DS line cut | fixed size = %g um', fixed_um(k)));
    if k == 1
        legend('Location', 'best');
    end
end
h = sgtitle('Extra: differential sedimentation line cuts');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig6);
saveas(fig6, fullfile(fig_dir, 'kernel_ds_linecuts_fixed_sizes.png'));

% -------- Extra: repo-law DS size-size compare --------
repo_ds_z = cell(numel(repo_law_list), 1);
repo_ds_all = [];
for k = 1:numel(repo_law_list)
    beta_day = ds_kernel_day(D1m, D2m, repo_law_list{k}, cfg, grid);
    repo_ds_z{k} = ds_plot_field(beta_day, D1m, D2m, 0.03);
    repo_ds_all = [repo_ds_all; repo_ds_z{k}(:)]; %#ok<AGROW>
end
repo_ds_cmax = max(abs(repo_ds_all));
repo_ds_levels = linspace(-repo_ds_cmax, repo_ds_cmax, 23);

fig7 = figure('Color', 'w', 'Position', [80 80 1500 1200]);
for k = 1:numel(repo_law_list)
    subplot(2, 2, k);
    contourf(D1m, D2m, repo_ds_z{k}, repo_ds_levels, 'LineColor', 'none');
    colormap(gca, div_cmap(256));
    clim([-repo_ds_cmax, repo_ds_cmax]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    xlabel('Particle 1 size (um)');
    ylabel('Particle 2 size (um)');
    title(sprintf('Extra DS | %s', repo_law_labels{k}));
    cb = colorbar;
    cb.Label.String = 'log_{10}(\beta_{DS}) [cm^3 day^{-1}]';
end
h = sgtitle('Extra: differential sedimentation kernel for repo settling laws');
set(h, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
plain_style(fig7);
saveas(fig7, fullfile(fig_dir, 'kernel_ds_size_size_repo_laws.png'));

copy_one(fullfile(fig_dir, 'kernel_brownian_temp_size_contours.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'kernel_shear_eps_size_contours.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'kernel_ds_size_size_settling_laws.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'kernel_ds_settling_speed_laws.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'kernel_ds_ratio_to_current.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'kernel_ds_linecuts_fixed_sizes.png'), docs_fig_dir);
copy_one(fullfile(fig_dir, 'kernel_ds_size_size_repo_laws.png'), docs_fig_dir);

fprintf('\nSaved kernel figures to: %s\n', fig_dir);

function beta_day = brownian_kernel_day(D1_um, D2_um, T_C, cfg)
[r1, ~] = image_diam_um_to_radii(D1_um, cfg);
[r2, ~] = image_diam_um_to_radii(D2_um, cfg);
r_pair = [r1(:)'; r2(:)'];
shape = reshape(KernelLibrary.brownian(r_pair, [], []), size(D1_um));
mu = cfg.kvisc * cfg.rho_fl;
T_K = T_C + 273.15;
conBr = (2.0/3.0) * cfg.k * T_K ./ mu;
beta_day = shape .* conBr * cfg.day_to_sec;
end

function beta_day = shear_kernel_day(D1_um, D2_um, eps_mks, cfg, grid)
[r1, rcons1] = image_diam_um_to_radii(D1_um, cfg);
[r2, rcons2] = image_diam_um_to_radii(D2_um, cfg);
r_pair = [r1(:)'; r2(:)'];
rcons_pair = [rcons1(:)'; rcons2(:)'];
param = struct('r_to_rg', cfg.r_to_rg, 'setcon', grid.setcon);
shape = reshape(KernelLibrary.curvilinearShear(r_pair, rcons_pair, param), size(D1_um));
eps_cgs = eps_mks * 1e4;
gamma = sqrt(eps_cgs / cfg.kvisc);
beta_day = shape .* gamma * cfg.day_to_sec;
end

function beta_day = ds_kernel_day(D1_um, D2_um, law_name, cfg, grid)
[r1, ~] = image_diam_um_to_radii(D1_um, cfg);
[r2, ~] = image_diam_um_to_radii(D2_um, cfg);
w1 = settling_speed_cms(D1_um, law_name, cfg, grid);
w2 = settling_speed_cms(D2_um, law_name, cfg, grid);
r_small = min(r1, r2) * cfg.r_to_rg;
beta_day = 0.5 * pi * abs(w1 - w2) .* (r_small .^ 2) * cfg.day_to_sec;
end

function w_mday = settling_speed_mday(D_um, law_name, cfg, grid)
w_mday = settling_speed_cms(D_um, law_name, cfg, grid) * cfg.day_to_sec / 100;
end

function w_cms = settling_speed_cms(D_um, law_name, cfg, grid)
D_cm = D_um * 1e-4;
[r_img, rcons] = image_diam_um_to_radii(D_um, cfg);
law = lower(string(law_name));
switch law
    case "current"
        w_cms = grid.setcon .* (rcons .^ 3) ./ max(r_img, realmin);
    case "kriest_8"
        w_mday = 66 .* (D_cm .^ 0.62);
        w_cms = (w_mday * 100) / cfg.day_to_sec;
    case "kriest_9"
        w_mday = 132 .* (D_cm .^ 0.62);
        w_cms = (w_mday * 100) / cfg.day_to_sec;
    case "siegel_2025"
        D_mm = D_cm * 10;
        w_mday = 20.2 .* (D_mm .^ 0.67);
        w_cms = (w_mday * 100) / cfg.day_to_sec;
    case "white"
        delta_rho = 0.0014 * ones(size(D_cm)); % g cm^-3, LC 2015 Fig. 1
        w_cms = white_drag_case_speed_cms(D_cm, delta_rho, cfg);
    case "laurenceau_cornec"
        D_mm = D_cm * 10;
        a_lc = 0.03;
        d3_lc = 1.8;
        rho_sol = 0.3 * 0.8 + 0.2 * 2.0 + 0.5 * 1.06;
        solid_frac = a_lc .* (D_mm .^ (d3_lc - 3.0));
        solid_frac = min(max(solid_frac, 0), 1);
        delta_rho = solid_frac .* max(rho_sol - cfg.rho_fl, 0);
        w_cms = white_drag_case_speed_cms(D_cm, delta_rho, cfg);
    otherwise
        error('Unknown law: %s', law_name);
end
end

function w_cms = white_drag_case_speed_cms(D_cm, delta_rho, cfg)
mu = cfg.rho_fl * cfg.kvisc;
w_cms = (cfg.g .* delta_rho .* (D_cm .^ 2)) ./ max(18 .* mu, realmin);
w_cms = max(w_cms, 0);
re0 = w_cms .* D_cm ./ max(cfg.kvisc, realmin);
mask = re0 > 0.5;
if any(mask(:))
    w_iter = w_cms(mask);
    d_iter = D_cm(mask);
    dr_iter = delta_rho(mask);
    for it = 1:30
        re = max(w_iter .* d_iter ./ max(cfg.kvisc, realmin), realmin);
        cd = 24 ./ re + 6 ./ (1 + sqrt(re)) + 0.4;
        w_new = sqrt(max((4/3) .* cfg.g .* dr_iter .* d_iter ./ max(cd .* cfg.rho_fl, realmin), 0));
        w_iter = 0.5 * (w_iter + w_new);
    end
    re_end = w_iter .* d_iter ./ max(cfg.kvisc, realmin);
    use_iter = re_end > 0.5;
    w_keep = w_cms(mask);
    w_keep(use_iter) = w_iter(use_iter);
    w_cms(mask) = w_keep;
end
end

function [r_img, rcons] = image_diam_um_to_radii(D_um, cfg)
D_cm = D_um * 1e-4;
grid0 = cfg.derive();
amfrac = grid0.amfrac;
fr_dim = cfg.fr_dim;
r_img = D_cm ./ (2.0 * cfg.r_to_rg);
v = (r_img ./ amfrac) .^ fr_dim;
rcons = ((0.75 / pi) .* v) .^ (1.0 / 3.0);
end

function copy_one(src_file, dst_dir)
[~, name, ext] = fileparts(src_file);
copyfile(src_file, fullfile(dst_dir, [name, ext]));
end

function lev = full_levels(zplot, nlev)
z = zplot(isfinite(zplot));
if isempty(z)
    lev = linspace(-2, 2, nlev);
    return
end
z_lo = min(z);
z_hi = max(z);
if z_hi <= z_lo
    z_hi = z_lo + 1;
end
lev = linspace(z_lo, z_hi, nlev);
end

function lev = zero_levels(zplot, nlev)
z = zplot(isfinite(zplot));
if isempty(z)
    lev = linspace(-1, 1, nlev);
    return
end
cmax = max(abs(z));
if ~(isfinite(cmax) && cmax > 0)
    cmax = 1;
end
lev = linspace(-cmax, cmax, nlev);
end

function lev = auto_levels(zplot, nlev)
z = zplot(isfinite(zplot));
if isempty(z)
    lev = linspace(-1, 1, nlev);
    return
end
if min(z) < 0 && max(z) > 0
    lev = zero_levels(zplot, nlev);
else
    lev = full_levels(zplot, nlev);
end
end

function cmap = mono_cmap(n, rgb_hi)
if nargin < 1
    n = 256;
end
if nargin < 2
    rgb_hi = [0.01 0.10 0.55];
end
t = linspace(0, 1, n)';
rgb_lo = [0.98 0.99 1.00];
cmap = (1 - t) .* rgb_lo + t .* rgb_hi;
end

function cmap = div_cmap(n)
if nargin < 1
    n = 256;
end
n1 = floor(n/2);
n2 = n - n1;
neg = [linspace(0.02, 1.00, n1)', linspace(0.08, 1.00, n1)', linspace(0.65, 1.00, n1)'];
pos = [ones(n2,1), linspace(1.00, 0.05, n2)', linspace(1.00, 0.05, n2)'];
cmap = [neg; pos];
end

function zplot = ds_plot_field(field_in, D1m, D2m, tol_log)
zplot = log10(max(field_in, realmin));
diag_mask = abs(log10(D1m) - log10(D2m)) < tol_log;
zgood = zplot(isfinite(zplot) & zplot > -100);
if isempty(zgood)
    zfloor = -12;
else
    zfloor = min(zgood) - 0.3;
end
zplot(zplot < zfloor) = zfloor;
n = size(zplot, 1);
for i = 1:n
    band = find(diag_mask(i, :));
    if isempty(band)
        continue
    end
    j1 = max(min(band) - 1, 1);
    j2 = min(max(band) + 1, size(zplot, 2));
    vals = [zplot(i, j1), zplot(i, j2)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vals = zplot(i, ~diag_mask(i, :));
        vals = vals(isfinite(vals));
    end
    if isempty(vals)
        fillv = zplot(i, band(1));
    else
        fillv = mean(vals);
    end
    zplot(i, band) = fillv;
end
end

function plain_style(fig)
ax = findall(fig, 'Type', 'axes');
for k = 1:numel(ax)
    set(ax(k), 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10);
    set(ax(k).Title,  'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
    set(ax(k).XLabel, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
    set(ax(k).YLabel, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
    if isprop(ax(k), 'ZLabel')
        set(ax(k).ZLabel, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
    end
end

lg = findall(fig, 'Type', 'legend');
for k = 1:numel(lg)
    set(lg(k), 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
end

cb = findall(fig, 'Type', 'ColorBar');
for k = 1:numel(cb)
    set(cb(k), 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10);
    set(cb(k).Label, 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'tex');
end

txt = findall(fig, 'Type', 'text');
for k = 1:numel(txt)
    set(txt(k), 'FontName', 'Calibri', 'FontWeight', 'normal', 'FontAngle', 'normal', 'Interpreter', 'tex');
    if strcmp(get(txt(k), 'Tag'), 'suptitle')
        set(txt(k), 'FontSize', 10);
    end
end
end
