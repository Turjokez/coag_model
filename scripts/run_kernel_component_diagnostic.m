% run_kernel_component_diagnostic.m
% Diagnose relative kernel contributions and alpha scaling.

clear; close all; clc;
setup_paths

cfg = SimulationConfig();
cfg.n_sections = 35;
cfg.alpha = 1.0;

grid = cfg.derive();
assembler = BetaAssembler(cfg, grid);

% Compute unscaled components
b_brown = assembler.computeFor('KernelBrown');
b_shear = assembler.computeFor('KernelCurSh');
b_ds    = assembler.computeFor('KernelCurDS');

% Scale components the same way as combineAndScale
b_brown_s = assembler.scaleBetas(b_brown, grid.conBr * cfg.day_to_sec);
b_shear_s = assembler.scaleBetas(b_shear, cfg.gamma * cfg.day_to_sec);
b_ds_s    = assembler.scaleBetas(b_ds,    grid.setcon * cfg.day_to_sec);

mag_brown = betaMagnitude(b_brown_s);
mag_shear = betaMagnitude(b_shear_s);
mag_ds    = betaMagnitude(b_ds_s);

total = mag_brown + mag_shear + mag_ds;
fprintf('Kernel contribution magnitudes (scaled, alpha=1):\n');
fprintf('  Brownian: %.3e (%.2f%%)\n', mag_brown, 100*mag_brown/total);
fprintf('  Shear:    %.3e (%.2f%%)\n', mag_shear, 100*mag_shear/total);
fprintf('  DS:       %.3e (%.2f%%)\n', mag_ds, 100*mag_ds/total);

% Check alpha scaling
cfg1 = copy(cfg);
cfg1.alpha = 1.0;
grid1 = cfg1.derive();
assembler1 = BetaAssembler(cfg1, grid1);
betas1 = assembler1.combineAndScale(b_brown, b_shear, b_ds);
mag1 = betaMagnitude(betas1);

cfg2 = copy(cfg);
cfg2.alpha = 1e4;
grid2 = cfg2.derive();
assembler2 = BetaAssembler(cfg2, grid2);
betas2 = assembler2.combineAndScale(b_brown, b_shear, b_ds);
mag2 = betaMagnitude(betas2);

fprintf('\nAlpha scaling check:\n');
fprintf('  mag(alpha=1)   = %.3e\n', mag1);
fprintf('  mag(alpha=1e4) = %.3e\n', mag2);
fprintf('  ratio = %.2e (expected ~1e4)\n', mag2/max(mag1, eps));

% --------- helper ---------
function mag = betaMagnitude(betas)
    mag = sum(abs(betas.b1(:))) + ...
          sum(abs(betas.b2(:))) + ...
          sum(abs(betas.b3(:))) + ...
          sum(abs(betas.b4(:))) + ...
          sum(abs(betas.b5(:))) + ...
          sum(abs(betas.b25(:)));
end
