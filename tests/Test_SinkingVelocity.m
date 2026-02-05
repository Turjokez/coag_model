classdef Test_SinkingVelocity < matlab.unittest.TestCase

methods (Test)

    function sinkingLawsAreFiniteAndMonotonic(tc)
        cfg = SimulationConfig();
        grid = cfg.derive();

        laws = {"current","kriest_8","kriest_9","siegel_2025","kriest_8_capped","kriest_8_flat"};
        % broad order-of-magnitude checks at ~1 mm (m/day)
        ranges = struct();
        ranges.current     = [1e-6, 1e6];
        ranges.kriest_8    = [1e-2, 1e3];
        ranges.kriest_9    = [1e-2, 2e3];
        ranges.siegel_2025 = [1e-2, 1e3];
        ranges.kriest_8_capped = [1e-2, 1e3];
        ranges.kriest_8_flat = [1e-2, 1e3];

        for i = 1:numel(laws)
            cfg.sinking_law = laws{i};
            v_cms = SettlingVelocityService.velocityForSections(grid, cfg);

            tc.verifyFalse(any(~isfinite(v_cms)));
            tc.verifyGreaterThanOrEqual(min(v_cms), 0);

            % monotonic non-decreasing with size (allow tiny numerical noise)
            dv = diff(v_cms(:));
            tol = 1e-12 * max(1, max(abs(v_cms)));
            tc.verifyGreaterThanOrEqual(min(dv), -tol);

            % order-of-magnitude check near ~1 mm
            D_mm = grid.getVolumeDiameters() * 10;
            [~, idx] = min(abs(D_mm - 1.0));
            v_mday = (v_cms(idx) / 100) * cfg.day_to_sec;

            r = ranges.(laws{i});
            tc.verifyGreaterThanOrEqual(v_mday, r(1));
            tc.verifyLessThanOrEqual(v_mday, r(2));
        end
    end

end
end
