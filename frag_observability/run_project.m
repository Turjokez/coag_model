% run_project.m
% One-command paper workflow for frag observability.

close all;

frag_root = fileparts(mfilename('fullpath'));
addpath(fullfile(frag_root, 'experiments'));
addpath(fullfile(frag_root, 'diagnostics'));
rehash;

opts = struct();
opts.sinking_laws = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
opts.apply_uvp = true;
opts.uvp_sample_L = 1.0;

fprintf('\nRunning frag observability paper workflow...\n');
fprintf('Sinking laws: %s\n', strjoin(opts.sinking_laws, ', '));
fprintf('UVP-like sampling: %s | sample volume = %.3g L\n', ...
    mat2str(opts.apply_uvp), opts.uvp_sample_L);

out = generate_first_pass_paper_figures(opts);

fprintf('\nDone.\n');
fprintf('Figures: %s\n', out.figure_dir);
fprintf('Tables:  %s\n', fullfile(frag_root, 'tables'));
fprintf('Notes:   %s\n', fullfile(frag_root, 'notes'));
