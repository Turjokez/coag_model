% run_gamma_sweep.m
% Sweep shear rate (gamma) and compare coag/sink ratio + export fractions.

clear; close all; clc;
setup_paths

gammas = [0.03, 0.1, 0.3, 1.0, 3.0]; % s^-1

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
base.box_depth = 2000;   % keep depth fixed for this sweep

% output folders
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'output', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

ratio_all = [];
D_mm = [];

fS_all = nan(1, numel(gammas));
fM_all = nan(1, numel(gammas));
fL_all = nan(1, numel(gammas));

for i = 1:numel(gammas)
    cfg = copy(base);
    cfg.gamma = gammas(i);

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

    fprintf('gamma=%.2f s^-1: small=%.3f  med=%.3f  large=%.3f\n', ...
        gammas(i), fS_all(i), fM_all(i), fL_all(i));
end

%% Figure: coag/sink ratio
fig1 = figure; hold on;
for i = 1:numel(gammas)
    loglog(D_mm, ratio_all(:,i), 'DisplayName', sprintf('\\gamma=%.2f', gammas(i)));
end
grid on;
xlabel('Diameter (mm)');
ylabel('Coag/Sink ratio');
title('Coagulation vs Sinking (day 50, \gamma sweep)');
legend('Location','best');
saveas(fig1, fullfile(fig_dir, 'gamma_coag_vs_sink.png'));

fprintf('\nSaved figure to: %s\n', fig_dir);
