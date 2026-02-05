classdef Test_Basics < matlab.unittest.TestCase

%i wil run by
% results = runtests("tests");
% table(results)

methods (Test)

    function kernelsAreFinite(tc)
        sim = CoagulationSimulation(SimulationConfig());
        B = buildBetas(sim);
        tc.verifyFalse(any(~isfinite(B.b1(:))));
        tc.verifyFalse(any(~isfinite(B.b25(:))));
    end

    function zeroStateGivesZeroRHS(tc)
        cfg = SimulationConfig();
        cfg.growth = 0; 
        cfg.c3     = 0;      % disagg loop off

        sim = CoagulationSimulation(cfg);

        B   = buildBetas(sim);
        L0  = zeros(cfg.n_sections);      % no growth/sinking in this test

        rhs = CoagulationRHS(B, L0, sim.operators.disagg_minus, sim.operators.disagg_plus, cfg);
        rhs.validate();

        dv  = rhs.evaluate(0, zeros(cfg.n_sections,1));

        tc.verifyFalse(any(isnan(dv)));
        tc.verifyLessThan(max(abs(dv)), 1e-14);
    end

    function coagulationConservesBiovolumeProxy(tc)
        cfg = SimulationConfig();
        cfg.growth = 0;
        cfg.c3     = 0;

        sim = CoagulationSimulation(cfg);

        B   = buildBetas(sim);
        L0  = zeros(cfg.n_sections);      % isolate coag

        rhs = CoagulationRHS(B, L0, sim.operators.disagg_minus, sim.operators.disagg_plus, cfg);
        rhs.validate();

        v0  = InitialSpectrumBuilder.initialSpectrum(cfg, sim.grid);
        dv  = rhs.evaluate(0, v0);

        leak = sum(dv .* sim.grid.av_vol(:)); % should be ~0 for coag-only
        tc.verifyFalse(isnan(leak));
        tc.verifyLessThan(abs(leak), 1e-10);
    end

    function sinkingRemovesInventoryProxy(tc)
        cfg = SimulationConfig();
        cfg.growth = 0;
        cfg.c3     = 0;

        sim = CoagulationSimulation(cfg);

        % Build sinking-only linear operator
        L = LinearProcessBuilder.linearMatrix(cfg, sim.grid);

        % Turn off coag by forcing betas = 0
        B = buildBetas(sim);
        B = zeroBetas(B);

        rhs = CoagulationRHS(B, L, sim.operators.disagg_minus, sim.operators.disagg_plus, cfg);
        rhs.validate();

        v0  = InitialSpectrumBuilder.initialSpectrum(cfg, sim.grid);
        dv  = rhs.evaluate(0, v0);

        dM = sum(dv .* sim.grid.av_vol(:));   % should be negative if sinking removes inventory
        tc.verifyFalse(isnan(dM));
        tc.verifyLessThan(dM, 0);
    end

end
end

% -------- helpers (kept outside class, small) ----------------------------
function B = buildBetas(sim)
b1 = sim.assembler.computeFor('KernelBrown');
b2 = sim.assembler.computeFor('KernelCurSh');
b3 = sim.assembler.computeFor('KernelCurDS');
B  = sim.assembler.combineAndScale(b1,b2,b3);
end

function B = zeroBetas(B)
B.b1(:)=0; B.b2(:)=0; B.b3(:)=0; B.b4(:)=0; B.b5(:)=0; B.b25(:)=0;
end