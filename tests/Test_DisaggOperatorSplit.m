classdef Test_DisaggOperatorSplit < matlab.unittest.TestCase

methods (Test)

    function redistributionMatchesResetLogic(tc)
        cfg = SimulationConfig();
        cfg.n_sections = 6;
        cfg.enable_disagg = true;
        cfg.disagg_mode = 'operator_split';
        cfg.disagg_frac_next = 2/3;
        cfg.disagg_C = 3.0;
        cfg.disagg_gamma = 0.15;

        grid = cfg.derive();
        d_low = 2.0 * (grid.amfrac * grid.v_lower.^grid.bmfrac);

        k = 4;
        d_max_cm = d_low(k);
        epsilon = (0.1 * cfg.disagg_C / d_max_cm)^(1/cfg.disagg_gamma);
        cfg.disagg_epsilon = epsilon;

        v = (1:6)';
        v_out = DisaggregationOperatorSplit.apply(v, grid, cfg, 0);

        expected = [1+2.5; 2+2.5; 3+10; 0; 0; 0];
        tc.verifyEqual(v_out, expected, 'AbsTol', 1e-12);
    end

    function dmaxAboveGridDoesNothing(tc)
        cfg = SimulationConfig();
        cfg.n_sections = 10;
        cfg.enable_disagg = true;
        cfg.disagg_mode = 'operator_split';
        cfg.disagg_dmax_cm = 1e6;

        grid = cfg.derive();
        v = (1:cfg.n_sections)';
        v_out = DisaggregationOperatorSplit.apply(v, grid, cfg, 0);

        tc.verifyEqual(v_out, v);
    end

end
end
