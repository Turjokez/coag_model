% run_slams_compare.m
% Compare SLAMS super-particle aggregation with sectional ODE model.
% Both models: 0-D slab, coagulation only, no sinking, no fragmentation.
% Same initial conditions, same kernels, same parameters.
%
% Key difference:
%   Sectional: solves Smoluchowski ODE via ode45 (continuous in time).
%   SLAMS-style: forward-Euler sub-steps on discrete super-particles;
%                merged aggregates snap to nearest Nn bin.
%
% State variable note:
%   The sectional model state is phi_k = N_k * av_vol(k) [volume concentration].
%   SLAMS works with N_k [#/cm^3] directly.
%   Conversion: N_k = phi_k / av_vol(k).

clear; clc;
this_file = mfilename('fullpath');
this_dir  = fileparts(this_file);
repo_root = fileparts(this_dir);
addpath(fullfile(repo_root, 'src'));

% --- Configuration ---
cfg = SimulationConfig( ...
    'n_sections',     20,           ...
    'sinking_law',    'kriest_8',   ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_sinking', false,        ...  % pure coagulation, no sinking
    'enable_disagg',  false,        ...
    'enable_coag',    true,         ...
    't_init',   0,                  ...
    't_final',  30,                 ...
    'delta_t',  1,                  ...
    'alpha',    1.0,                ...
    'gamma',    0.1,                ...
    'num_1',    1e3                 ...
);
cfg.validate();
grid = cfg.derive();

% Initial spectrum from sectional model builder
% phi0(k) = N_k * av_vol(k)  [volume concentration, cm^3/cm^3]
phi0 = InitialSpectrumBuilder.initialSpectrum(cfg, grid);

% SLAMS bin mass basis: Nn = 2^(k-1), volume per aggregate = Nn * v0
Nn_bins = 2 .^ (0:cfg.n_sections-1)';
Vagg_sl = Nn_bins * grid.v0;             % [cm^3 per aggregate]

% Convert to number concentration for SLAMS [#/cm^3]
n0_num = phi0 ./ Vagg_sl;

% --- Run sectional model ---
fprintf('Running sectional model...\n');
sim    = CoagulationSimulation(cfg);
result = sim.run('v0', phi0);

t_sec   = result.time;                   % [days]
phi_sec = result.concentrations;         % [n_times x n_sec], volume concentration
N_sec   = phi_sec ./ Vagg_sl';           % [n_times x n_sec, #/cm^3]
r_sec   = grid.getFractalRadii() * 1e4;  % [um]

% --- Run SLAMS-style model ---
fprintf('Running SLAMS-style model...\n');
[t_sl, Nn_arr, N_sl] = slams_run(cfg, grid, n0_num);

a0    = grid.a0;
D_sl  = 2.0;
r_sl  = a0 * Nn_arr .^ (1/D_sl) * 1e4;  % [um]

% --- Common x-axis: volume-equivalent diameter [um] ---
% Sectional: use conserved-volume radius (fr_dim=2.33)
% SLAMS:     d_eq = 2 * a0 * Nn^(1/3)  (D=2)
d_sec = 2 * grid.getConservedRadii() * 1e4;   % [um]
d_sl  = 2 * a0 * Nn_arr .^ (1/3) * 1e4;       % [um]

% Shared axis limits
plot_days = [0, 5, 10, 20, 30];
colors    = lines(length(plot_days));

all_N = [N_sec(N_sec > 0); N_sl(N_sl > 0)];
ylim_vals = [min(all_N)*0.5, max(all_N)*2];
xlim_vals = [min([d_sec; d_sl])*0.8, max([d_sec; d_sl])*1.2];

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% --- Figure 1: N(d_eq) at selected times ---
f1 = figure;
subplot(1,2,1);
hold on;
for k = 1:length(plot_days)
    td = plot_days(k);
    it = find(abs(t_sec - td) == min(abs(t_sec - td)), 1);
    semilogy(d_sec, N_sec(it,:), '-', 'Color', colors(k,:), ...
             'DisplayName', sprintf('t=%d d', td));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(xlim_vals); ylim(ylim_vals);
xlabel('d_{eq} [\mum]');
ylabel('N [#/cm^3]');
title('Sectional ODE');
legend show;

subplot(1,2,2);
hold on;
for k = 1:length(plot_days)
    td = plot_days(k);
    it = find(abs(t_sl - td) == min(abs(t_sl - td)), 1);
    semilogy(d_sl, N_sl(it,:), '-', 'Color', colors(k,:), ...
             'DisplayName', sprintf('t=%d d', td));
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim(xlim_vals); ylim(ylim_vals);
xlabel('d_{eq} [\mum]');
ylabel('N [#/cm^3]');
title('SLAMS-style (D=2)');
legend show;

sgtitle('N(d_{eq}): sectional vs SLAMS-style (0-D, coag only)');
saveas(f1, fullfile(fig_dir, 'may13_0d_spectra.png'));

% --- Figure 2: total N vs time ---
f2 = figure;
semilogy(t_sec, sum(N_sec, 2), 'b-',  'DisplayName', 'Sectional');
hold on;
semilogy(t_sl,  sum(N_sl,  2), 'r--', 'DisplayName', 'SLAMS');
hold off;
xlabel('time [days]');
ylabel('total N [#/cm^3]');
legend show;
title('Total particle count vs time');
saveas(f2, fullfile(fig_dir, 'may13_0d_totalN.png'));

% --- Figure 3: volume spectrum at t=0 and t=30, both models ---
% phi_k = N_k * v_agg(k) [cm^3/cm^3] — where mass actually is
it0_sec  = 1;
it30_sec = find(abs(t_sec - 30) == min(abs(t_sec - 30)), 1);
it30_sl  = find(abs(t_sl  - 30) == min(abs(t_sl  - 30)), 1);
v_primary = (4/3) * pi * a0^3;             % volume of one primary particle [cm^3]
phi_sl_30 = N_sl(it30_sl, :)' .* (Nn_arr * v_primary); % SLAMS volume spectrum, t=30

% Mask near-zero SLAMS bins — prevents y-axis collapsing to 100 decades
phi_sl_30_plot = phi_sl_30;
phi_sl_30_plot(phi_sl_30_plot < 1e-30) = NaN;

f3 = figure;
loglog(d_sec, phi_sec(it0_sec, :),  'k-',  'DisplayName', 't=0 (both)');
hold on;
loglog(d_sec, phi_sec(it30_sec, :), 'b-',  'DisplayName', 'Sectional t=30');
loglog(d_sl,  phi_sl_30_plot,       'r--', 'DisplayName', 'SLAMS t=30');
hold off;
xlabel('d_{eq} [\mum]');
ylabel('\phi [cm^3/cm^3]');
legend show;
title('Volume spectrum: sectional vs SLAMS (0-D, t=30 d)');
saveas(f3, fullfile(fig_dir, 'may13_0d_volspec.png'));

% --- Volume conservation check ---
% total particle volume should be conserved in pure coagulation.
% For sectional: total volume = sum(phi_k).
% For SLAMS: total volume = sum(N_k * Nn_k * v_primary).
V_agg     = Nn_arr * v_primary;             % volume per aggregate [cm^3]

vol_sec = sum(phi_sec, 2);                  % total volume [cm^3/cm^3]
vol_sl  = N_sl * V_agg;                     % total volume [cm^3/cm^3]

fprintf('\nVolume conservation check (total particle volume):\n');
fprintf('  Sectional: initial=%.3e, final=%.3e, change=%.2f%%\n', ...
    vol_sec(1), vol_sec(end), 100*(vol_sec(end)-vol_sec(1))/vol_sec(1));
fprintf('  SLAMS:     initial=%.3e, final=%.3e, change=%.2f%%\n', ...
    vol_sl(1),  vol_sl(end),  100*(vol_sl(end) -vol_sl(1)) /vol_sl(1));

fprintf('\nNumber check (total particle count, should decrease with coagulation):\n');
fprintf('  Sectional: initial=%.3e, final=%.3e, change=%.2f%%\n', ...
    sum(N_sec(1,:)), sum(N_sec(end,:)), ...
    100*(sum(N_sec(end,:))-sum(N_sec(1,:)))/sum(N_sec(1,:)));
fprintf('  SLAMS:     initial=%.3e, final=%.3e, change=%.2f%%\n', ...
    sum(N_sl(1,:)), sum(N_sl(end,:)), ...
    100*(sum(N_sl(end,:))-sum(N_sl(1,:)))/sum(N_sl(1,:)));
