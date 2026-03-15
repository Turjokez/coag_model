# How To Run

This paper workspace is meant to stay inside `frag_observability/`.

It does not change the core model physics.
It reads the current box-model code from `coag_model_final/`.

## Fastest way

From MATLAB, at the repo root, use the root-level wrapper:

```matlab
run('run_frag_observability_project.m')
```

Or, if you want to call the paper runner directly, only do this when your current folder is the repo root:

```matlab
run('frag_observability/run_project.m')
```

If your current folder is already `frag_observability/`, use:

```matlab
run('run_project.m')
```

That runs the full paper-side workflow:

- experiment matrix
- summary tables
- standard diagnostic figures
- first-pass paper figures
- note update files

## What the default automatic run does

The default automatic run uses:

- 4 sinking laws:
  - `current`
  - `kriest_8`
  - `kriest_9`
  - `siegel_2025`
- 4 fragmentation levels:
  - `no_frag`
  - `weak`
  - `medium`
  - `strong`
- image and volume diameter outputs
- UVP-like Poisson sampling:
  - `apply_uvp = true`
  - `uvp_sample_L = 1.0`

Default model count:

- for each sinking law:
  - 3 baseline scale-scan runs
  - 4 main paper cases
- total default model runs:
  - `4 * (3 + 4) = 28`

## Main outputs

Tables are written to:

- `frag_observability/tables/frag_observability_case_summary.csv`
- `frag_observability/tables/frag_observability_baseline_scales.csv`
- `frag_observability/tables/frag_observability_case_summary.mat`

Paper figures are written to:

- `frag_observability/figures/fig01_b_kappa_separability.png`
- `frag_observability/figures/fig02_baseline_psd_fit_residuals.png`
- `frag_observability/figures/fig03_kappa_vs_fragmentation_strength.png`
- `frag_observability/figures/fig04_clean_window_noisy_psd.png`
- `frag_observability/figures/fig05_delta_b_delta_kappa_relative.png`
- `frag_observability/figures/fig06_noise_1L_vs_10L_delta_space.png`

Framework figures are also written to the same folder.

Notes are written to:

- `frag_observability/notes/first_pass_figure_note.md`
- `frag_observability/notes/initial_findings.md`

## If you want only the experiment tables and standard framework figures

From MATLAB:

```matlab
addpath('frag_observability/experiments');
opts = struct();
opts.sinking_laws = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
opts.apply_uvp = false;
run_frag_observability_framework(opts);
```

## If you want the paper figures directly

From MATLAB:

```matlab
addpath('frag_observability/experiments');
opts = struct();
opts.sinking_laws = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
opts.apply_uvp = true;
opts.uvp_sample_L = 1.0;
generate_first_pass_paper_figures(opts);
```

## Command-line automatic run

This machine does **not** have `matlab` on `PATH`, but this direct app path works from the repo root:

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "cd('/Users/turjo/Documents/coag_model'); run('run_frag_observability_project.m')"
```

If `matlab` is on your shell path, from the repo root, this shorter form also works:

```bash
matlab -batch "run('run_frag_observability_project.m')"
```

If you already have saved tables and only want to rerender the paper figures:

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "cd('/Users/turjo/Documents/coag_model'); addpath('frag_observability/experiments'); generate_first_pass_paper_figures(struct('use_saved_results', true));"
```

## Notes

- The workflow saves outputs only inside `frag_observability/`.
- Re-running will overwrite paper-side files with the same names in this folder.
- Core files in `coag_model_final/src/` are not edited by this run.
- The MATLAB batch run may print a warning about a missing old startup folder. The workflow still runs fine after that warning.
