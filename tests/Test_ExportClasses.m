% Test_ExportClasses.m


classdef Test_ExportClasses < matlab.unittest.TestCase

methods (Test)

    function exportVectorExistsAndFinite(tc)
        cfg = SimulationConfig();
        cfg.growth = 0;
        cfg.c3     = 0;

        sim = CoagulationSimulation(cfg);
        r   = sim.run('tspan', 0:1:2);

        tc.verifyTrue(isfield(r, 'output_data'));
        od = r.output_data;

        % export flux must exist
        tc.verifyTrue(isfield(od, 'fluxsect'));

        % model-consistent diameter (cm)
        d_cm = sim.grid.getVolumeDiameters();
        d_cm = d_cm(:);

        flux = od.fluxsect(end, :).';

        tc.verifyEqual(numel(d_cm), numel(flux));
        tc.verifyFalse(any(~isfinite(d_cm)));
        tc.verifyFalse(any(~isfinite(flux)));
    end

    function exportFractionsSumToOneWhenExportPositive(tc)
        cfg = SimulationConfig();
        cfg.growth = 0;
        cfg.c3     = 0;

        sim = CoagulationSimulation(cfg);
        r   = sim.run('tspan', 0:1:2);

        od = r.output_data;

        % model-consistent diameter (cm)
        d_cm = sim.grid.getVolumeDiameters();
        d_cm = d_cm(:);

        flux = od.fluxsect(end, :).';

        [fS,fM,fL,tot] = export_size_classes(d_cm, flux);

        if tot.total > 0
            tc.verifyLessThan(abs((fS+fM+fL) - 1), 1e-10);
            tc.verifyGreaterThanOrEqual(fS, 0);
            tc.verifyGreaterThanOrEqual(fM, 0);
            tc.verifyGreaterThanOrEqual(fL, 0);
        else
            tc.verifyEqual(fS+fM+fL, 0);
        end
    end

end
end