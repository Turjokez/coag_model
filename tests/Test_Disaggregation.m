classdef Test_Disaggregation < matlab.unittest.TestCase

methods (Test)

    function zeroInputReturnsZero(tc)
        cfg = SimulationConfig();
        v = zeros(cfg.n_sections, 1);
        term = Disaggregation.netTerm(v, cfg);
        tc.verifyEqual(term, zeros(cfg.n_sections,1));
    end

    function negativeInputIsClipped(tc)
        cfg = SimulationConfig();
        cfg.c3 = 0.2; cfg.c4 = 1.45;

        v = -abs(linspace(1, 0.1, cfg.n_sections))';
        term_neg  = Disaggregation.netTerm(v, cfg);
        term_clip = Disaggregation.netTerm(max(v, 0), cfg);

        tc.verifyEqual(term_neg, term_clip);
        tc.verifyFalse(any(isnan(term_neg)));
        tc.verifyFalse(any(isinf(term_neg)));
    end

    function positiveInputIsNonzero(tc)
        cfg = SimulationConfig();
        cfg.c3 = 0.2; cfg.c4 = 1.45;

        v = linspace(1, 0.1, cfg.n_sections)';
        term = Disaggregation.netTerm(v, cfg);

        tc.verifyGreaterThan(max(abs(term)), 0);
    end

end
end

