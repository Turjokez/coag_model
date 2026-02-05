% setup_paths.m  (SCRIPT, not a function)
% Adds project folders to the MATLAB path safely (no recursion).

% Project root = folder where this setup_paths.m lives
root = fileparts(mfilename('fullpath'));

% Add only the main folders you need
addpath(root);
addpath(fullfile(root,'src'));
addpath(fullfile(root,'scripts'));
addpath(fullfile(root,'diagnostics'));
addpath(fullfile(root,'data'));

% Refresh MATLAB path cache
rehash;

% Optional: show confirmation
fprintf('setup_paths: added paths under %s\n', root);