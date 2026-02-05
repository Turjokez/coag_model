clear; close all;
setup_paths

cfg = SimulationConfig();
sim = CoagulationSimulation(cfg);
res = sim.run();

d_cm_v = res.output_data.diam_v(:);             % volume diameter [cm]
d_cm_i = res.output_data.diam_i(:);             % image diameter [cm]
Fv_final = res.output_data.fluxspec(end,:).';   % volume flux spectrum
Fi_final = res.output_data.fluxspec_i(end,:).'; % image flux spectrum

[fS_v,fM_v,fL_v,tot_v] = export_size_classes(d_cm_v, Fv_final, 'volume');
[fS_i,fM_i,fL_i,tot_i] = export_size_classes(d_cm_i, Fi_final, 'image');

fprintf('\n=== Export fractions at final time ===\n');
fprintf('Volume diam: small=%.4f  med=%.4f  large=%.4f  (Total=%.6g)\n', ...
    fS_v, fM_v, fL_v, tot_v.total);
fprintf('Image  diam: small=%.4f  med=%.4f  large=%.4f  (Total=%.6g)\n', ...
    fS_i, fM_i, fL_i, tot_i.total);

figure;
bar([fS_v fM_v fL_v]);
set(gca,'XTickLabel',{'Small','Medium','Large'});
ylabel('Fraction of total export');
title('Export size-class fractions (volume diam)');
grid on;

figure;
bar([fS_i fM_i fL_i]);
set(gca,'XTickLabel',{'Small','Medium','Large'});
ylabel('Fraction of total export');
title('Export size-class fractions (image diam)');
grid on;

figure;
F = max(Fv_final,0);
loglog(d_cm_v*10, F, 'o-');
xline(0.5,'--'); xline(2.0,'--');
xlabel('Diameter (mm)'); ylabel('Export flux spectrum');
title('Final-time export by size (volume diam)');
grid on;
