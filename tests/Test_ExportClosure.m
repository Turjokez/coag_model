classdef Test_ExportClosure < matlab.unittest.TestCase

methods (Test)

    function sinkingLossMatchesRHS(tc)
        cfg = SimulationConfig();
        cfg.enable_coag   = false;
        cfg.enable_disagg = false;
        cfg.enable_linear = true;
        cfg.enable_sinking = true;
        cfg.growth = 0;

        sim = CoagulationSimulation(cfg);

        % Build sinking-only linear operator
        L = LinearProcessBuilder.linearMatrix(cfg, sim.grid);

        % Turn off coag by forcing betas = 0
        B = buildZeroBetas(cfg.n_sections);

        rhs = CoagulationRHS(B, L, sim.operators.disagg_minus, sim.operators.disagg_plus, cfg);
        rhs.validate();

        v0  = InitialSpectrumBuilder.initialSpectrum(cfg, sim.grid);
        dv  = rhs.evaluate(0, v0);

        % OutputGenerator sink loss (state/day) should match RHS sinking term
        t = 0;
        Y = v0(:)';
        od = OutputGenerator.spectraAndFluxes(t, Y, sim.grid, cfg);
        sinkLoss = od.sinkLossSect(1, :).';

        diff = dv + sinkLoss; % RHS should be -sinkLoss
        tol = 1e-10 * max(1, max(abs(dv)));

        tc.verifyLessThanOrEqual(max(abs(diff)), tol);
        tc.verifyLessThan(sum(dv .* sim.grid.av_vol(:)), 0);
    end

end
end

% -------- helpers ----------------------------
function B = buildZeroBetas(n)
B = BetaMatrices();
B.b1  = zeros(n,n);
B.b2  = zeros(n,n);
B.b3  = zeros(n,n);
B.b4  = zeros(n,n);
B.b5  = zeros(n,n);
B.b25 = zeros(n,n);
end
