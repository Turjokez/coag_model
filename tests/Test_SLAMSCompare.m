classdef Test_SLAMSCompare < matlab.unittest.TestCase
% Small tests for 0-D SLAMS vs sectional comparison.

methods (Test)

    function slamsVolumeConservation(tc)
        [cfg, grid, n0_num, Vagg] = setupCase();
        [~, ~, N_sl] = slams_run(cfg, grid, n0_num);

        vol0 = sum(N_sl(1,:)   .* Vagg');
        volf = sum(N_sl(end,:) .* Vagg');
        rel  = abs((volf - vol0) / max(vol0, eps));

        tc.verifyLessThan(rel, 1e-12);
    end

    function sectionalTopBinTruncationSignature(tc)
        [loss20, top0_20, topf_20] = runSectional(20);
        [loss40, ~, ~] = runSectional(40);

        % Initial spectrum is not concentrated near upper bound.
        tc.verifyLessThan(top0_20, 1e-10);
        % During run, material reaches top bin.
        tc.verifyGreaterThan(topf_20, 1e-3);
        % More sections reduce truncation loss.
        tc.verifyLessThan(loss20, loss40);
    end

    function dAndSnapBothChangeResult(tc)
        [cfg, grid, n0_num, ~] = setupCase();

        optsA = struct('D', 2.00, 'n_sub', 24, 'snap_mode', 'nearest');
        optsB = struct('D', 2.33, 'n_sub', 24, 'snap_mode', 'nearest');
        optsC = struct('D', 2.33, 'n_sub', 24, 'snap_mode', 'split');

        [~, ~, NA] = slams_run(cfg, grid, n0_num, optsA);
        [~, ~, NB] = slams_run(cfg, grid, n0_num, optsB);
        [~, ~, NC] = slams_run(cfg, grid, n0_num, optsC);

        NendA = sum(NA(end,:));
        NendB = sum(NB(end,:));
        NendC = sum(NC(end,:));

        % D change effect at fixed snap mode.
        tc.verifyGreaterThan(abs(NendB - NendA) / max(NendA, eps), 0.01);
        % Snap-mode effect at fixed D.
        tc.verifyGreaterThan(abs(NendC - NendB) / max(NendB, eps), 0.005);
    end

    function slamsSubstepConvergence24to48(tc)
        [cfg, grid, n0_num, ~] = setupCase();
        opts24 = struct('D', 2.00, 'n_sub', 24, 'snap_mode', 'nearest');
        opts48 = struct('D', 2.00, 'n_sub', 48, 'snap_mode', 'nearest');

        [~, ~, N24] = slams_run(cfg, grid, n0_num, opts24);
        [~, ~, N48] = slams_run(cfg, grid, n0_num, opts48);

        Nend24 = sum(N24(end,:));
        Nend48 = sum(N48(end,:));
        rel = abs(Nend24 - Nend48) / max(Nend48, eps);

        tc.verifyLessThan(rel, 0.01);
    end

end
end

function [cfg, grid, n0_num, Vagg] = setupCase()
this_test = mfilename('fullpath');
tests_dir = fileparts(this_test);
repo_root = fileparts(tests_dir);
addpath(fullfile(repo_root, 'slams_compare'));

cfg = SimulationConfig( ...
    'n_sections', 20, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_sinking', false, ...
    'enable_disagg', false, ...
    'enable_coag', true, ...
    't_init', 0, ...
    't_final', 30, ...
    'delta_t', 1, ...
    'alpha', 1.0, ...
    'gamma', 0.1, ...
    'num_1', 1e3);
grid = cfg.derive();
phi0 = InitialSpectrumBuilder.initialSpectrum(cfg, grid);
Nn   = 2 .^ (0:cfg.n_sections-1)';
Vagg = Nn * grid.v0;
n0_num = phi0 ./ Vagg;
end

function [volChangePct, topFracInit, topFracFinal] = runSectional(nsec)
cfg = SimulationConfig( ...
    'n_sections', nsec, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_sinking', false, ...
    'enable_disagg', false, ...
    'enable_coag', true, ...
    't_init', 0, ...
    't_final', 30, ...
    'delta_t', 1, ...
    'alpha', 1.0, ...
    'gamma', 0.1, ...
    'num_1', 1e3);
grid = cfg.derive();
phi0 = InitialSpectrumBuilder.initialSpectrum(cfg, grid);
sim = CoagulationSimulation(cfg);
res = sim.run('v0', phi0);
phiF = res.concentrations(end,:)';

volChangePct = 100 * (sum(phiF) - sum(phi0)) / sum(phi0);
topFracInit = phi0(end) / sum(phi0);
topFracFinal = phiF(end) / max(sum(phiF), eps);
end
