% run_compare_current_vs_kriest8.m
% One matched 0-D case. Only sinking law changes.

clear; close all; clc;
setup_paths

laws = {'kriest_8','current'};

% Same setup as current main diagnostics
base = SimulationConfig();
base.t_final = 200;
base.delta_t = 2;
base.n_sections = 35;
base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;
base.alpha = 1.0;
base.enable_coag = true;
base.enable_linear = true;
base.enable_sinking = true;
base.enable_disagg = true;
base.c3 = 0.02;
base.r_to_rg = 1.6;
base.growth = 0;
base.dz = 65;
base.box_depth = 2000;
base.sinking_size = 'volume';
base.sinking_scale = 1.0;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
tab_dir = fullfile(project_root, 'output', 'tables');
log_dir = fullfile(project_root, 'output', 'logs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
if ~exist(tab_dir,'dir'), mkdir(tab_dir); end
if ~exist(log_dir,'dir'), mkdir(log_dir); end

t = [];
d_um = [];
set_vel = [];
total_flux = [];
total_mass = [];
final_nspec = [];
summary = struct([]);

for i = 1:numel(laws)
    cfg = copy(base);
    cfg.sinking_law = laws{i};

    sim = CoagulationSimulation(cfg);
    res = sim.run();
    if isempty(t)
        t = res.time(:);
        d_um = res.output_data.diam_v(:) * 1e4;
        set_vel = zeros(numel(d_um), numel(laws));
        total_flux = zeros(numel(t), numel(laws));
        total_mass = zeros(numel(t), numel(laws));
        final_nspec = zeros(numel(d_um), numel(laws));
    end

    set_vel(:,i) = res.output_data.set_vel(:);
    total_flux(:,i) = res.output_data.total_flux(:);
    total_mass(:,i) = res.output_data.total_mass(:);
    final_nspec(:,i) = res.output_data.nspec_v(end,:).';

    [~, idx_1mm] = min(abs(d_um - 1000));
    [peak_flux, idx_peak] = max(res.output_data.total_flux(:));
    [fS, fM, fL] = export_size_classes(res.output_data.diam_v(:), ...
        res.output_data.fluxspec(end,:).', 'volume');

    summary(i).law = laws{i};
    summary(i).sink_1mm_mday = res.output_data.set_vel(idx_1mm);
    summary(i).sink_max_mday = max(res.output_data.set_vel);
    summary(i).final_inventory = res.output_data.total_mass(end);
    summary(i).final_total_flux = res.output_data.total_flux(end);
    summary(i).peak_flux = peak_flux;
    summary(i).peak_flux_day = t(idx_peak);
    summary(i).small_frac = fS;
    summary(i).medium_frac = fM;
    summary(i).large_frac = fL;
    summary(i).min_state = min(res.concentrations(:));
end

fig = figure('Color','w','Position',[100 100 980 720]);

subplot(2,2,1)
loglog(d_um, set_vel(:,1), 'k', 'LineWidth', 1.5); hold on
loglog(d_um, set_vel(:,2), 'r', 'LineWidth', 1.5);
xlabel('Diameter (um)');
ylabel('Sinking speed (m day^{-1})');
legend(laws, 'Location', 'northwest');
title('Sinking speed');

subplot(2,2,2)
plot(t, total_flux(:,1), 'k', 'LineWidth', 1.5); hold on
plot(t, total_flux(:,2), 'r', 'LineWidth', 1.5);
xlabel('Time (day)');
ylabel('Total export flux');
legend(laws, 'Location', 'best');
title('Export flux');

subplot(2,2,3)
plot(t, total_mass(:,1), 'k', 'LineWidth', 1.5); hold on
plot(t, total_mass(:,2), 'r', 'LineWidth', 1.5);
xlabel('Time (day)');
ylabel('Inventory in box');
legend(laws, 'Location', 'best');
title('Inventory');

subplot(2,2,4)
loglog(d_um, max(final_nspec(:,1), eps), 'k', 'LineWidth', 1.5); hold on
loglog(d_um, max(final_nspec(:,2), eps), 'r', 'LineWidth', 1.5);
xlabel('Diameter (um)');
ylabel('Final number spectrum');
legend(laws, 'Location', 'best');
title('Final spectrum');

exportgraphics(fig, fullfile(fig_dir, 'current_vs_kriest8_matched_case.png'), ...
    'Resolution', 200);

T = table( ...
    string({summary.law})', ...
    [summary.sink_1mm_mday]', ...
    [summary.sink_max_mday]', ...
    [summary.final_inventory]', ...
    [summary.final_total_flux]', ...
    [summary.peak_flux]', ...
    [summary.peak_flux_day]', ...
    [summary.small_frac]', ...
    [summary.medium_frac]', ...
    [summary.large_frac]', ...
    [summary.min_state]', ...
    'VariableNames', { ...
    'law','sink_1mm_mday','sink_max_mday','final_inventory', ...
    'final_total_flux','peak_flux','peak_flux_day', ...
    'small_frac','medium_frac','large_frac','min_state'});

writetable(T, fullfile(tab_dir, 'current_vs_kriest8_matched_case.csv'));

fid = fopen(fullfile(log_dir, 'current_vs_kriest8_matched_case.txt'), 'w');
fprintf(fid, 'Matched case: only sinking law changes\n');
fprintf(fid, 't_final=%g day, delta_t=%g day, n_sections=%g\n', ...
    base.t_final, base.delta_t, base.n_sections);
fprintf(fid, 'PP on, coag on, sinking on, disagg on, c3=%g, r_to_rg=%g, box_depth=%g m\n\n', ...
    base.c3, base.r_to_rg, base.box_depth);
for i = 1:height(T)
    fprintf(fid, '%s\n', T.law{i});
    fprintf(fid, '  sink at 1 mm      : %.6g m/day\n', T.sink_1mm_mday(i));
    fprintf(fid, '  max sink speed    : %.6g m/day\n', T.sink_max_mday(i));
    fprintf(fid, '  final inventory   : %.6g\n', T.final_inventory(i));
    fprintf(fid, '  final total flux  : %.6g\n', T.final_total_flux(i));
    fprintf(fid, '  peak flux         : %.6g at day %.6g\n', T.peak_flux(i), T.peak_flux_day(i));
    fprintf(fid, '  export frac S/M/L : %.6f  %.6f  %.6f\n', ...
        T.small_frac(i), T.medium_frac(i), T.large_frac(i));
    fprintf(fid, '  min state         : %.6g\n', T.min_state(i));
    fprintf(fid, '\n');
end
fclose(fid);

disp(T)
disp('Saved figure, table and log.')
