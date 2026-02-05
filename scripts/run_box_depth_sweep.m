% run_box_depth_sweep.m
% Sweep box_depth values and compare coag/sink ratio + export fractions.

clear; close all; clc;
setup_paths

depths = [65, 500, 2000, 5000]; % meters

% --------- base configuration ---------
base = SimulationConfig();
base.t_final = 50;
base.delta_t = 1;
base.n_sections = 35;

base.enable_pp = true;
base.pp_bin = 1;
base.pp_source = 1e-8;

base.enable_coag = true;
base.enable_sinking = true;
base.enable_disagg = false;
base.enable_linear = true;
base.growth = 0;

base.sinking_law = 'kriest_8';
base.sinking_size = 'volume';
base.sinking_scale = 1.0;

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

ratio_all = [];
D_mm = [];

fS_all = nan(1, numel(depths));
fM_all = nan(1, numel(depths));
fL_all = nan(1, numel(depths));

for i = 1:numel(depths)
    cfg = copy(base);
    cfg.box_depth = depths(i);

    sim = CoagulationSimulation(cfg);
    res = sim.run();

    if isempty(D_mm)
        gr = cfg.derive();
        D_mm = gr.getVolumeDiameters() * 10;
    end

    v = res.concentrations(end,:)';
    [term1, term2, term3, ~, ~] = sim.rhs.rateTerms(v);
    coag = abs(term1 + term2);
    sink = abs(term3);
    ratio_all(:,i) = coag ./ max(sink, eps);

    % export fractions at day 50
    F = res.output_data.fluxspec;
    if size(F,1) == numel(res.time)
        flux_final = F(end,:).';
    else
        flux_final = F(:,end);
    end
    [fS_all(i), fM_all(i), fL_all(i)] = export_size_classes(D_mm/10, flux_final); % mm -> cm

    fprintf('box_depth=%g m: small=%.3f  med=%.3f  large=%.3f\n', ...
        depths(i), fS_all(i), fM_all(i), fL_all(i));
end

%% Figure: coag/sink ratio
fig1 = figure; hold on;
for i = 1:numel(depths)
    loglog(D_mm, ratio_all(:,i), 'DisplayName', sprintf('H=%g m', depths(i)));
end
grid on;
xlabel('Diameter (mm)');
ylabel('Coag/Sink ratio');
title('Coagulation vs Sinking (day 50, box\_depth sweep)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'boxdepth_coag_vs_sink.png'));

fprintf('\nSaved figure to: %s\n', fig_dir);
