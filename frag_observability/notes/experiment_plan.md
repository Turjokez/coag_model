# Experiment Plan

This is the small paper experiment plan for frag observability.

## Goal

Use the current repo model to ask:

- can fragmentation be seen in PSD shape
- is `kappa` more useful than slope `b`
- can we still separate cases under a small core set of sinking laws

## Cases

The framework runs 4 frag cases:

- no fragmentation
- weak fragmentation
- medium fragmentation
- strong fragmentation

The default frag strengths use the same values already used in the repo:

- weak: `epsilon = 1e-8`
- medium: `epsilon = 1e-6`
- strong: `epsilon = 1e-4`

## Sinking laws

The framework now uses a 4-law core that already exists in the repo:

- `current`
- `kriest_8`
- `kriest_9`
- `siegel_2025`

Why these 4:

- they are part of the main repo sinking-law list
- they are already supported directly by `SimulationConfig` and `SettlingVelocityService`
- they give a broader comparison than a 2-law test without making the first paper pass too large

For each law, the framework first scans a small set of sinking scales and picks the no-frag case with the best image-space power-law fit.

## What is measured

For each run, the framework computes:

- slope `b` from `log10(N)` vs `log10(D)`
- power-law residuals
- quadratic residual curvature `kappa = c2`

It does this for:

- image diameter output
- volume diameter output

## Optional sampling

If `apply_uvp = true`, the framework also applies a simple Poisson sample step.

This is closest to the repo UVP-like logic in image space, but the code also allows the same sample step in volume space so both outputs stay parallel.

Default sample setting:

- `uvp_sample_L = 1.0`

## Output files

Tables go to:

- `frag_observability/tables/frag_observability_case_summary.csv`
- `frag_observability/tables/frag_observability_baseline_scales.csv`
- `frag_observability/tables/frag_observability_case_summary.mat`

Figures go to:

- `frag_observability/figures/frag_observability_psd_with_fit.png`
- `frag_observability/figures/frag_observability_residuals.png`
- `frag_observability/figures/frag_observability_b_kappa.png`
- `frag_observability/figures/frag_observability_b_kappa_sampled.png` if sampling is on

## Main runner

Run:

```matlab
run_frag_observability_framework
```

Or with sampling on:

```matlab
opts.apply_uvp = true;
opts.uvp_sample_L = 1.0;
run_frag_observability_framework(opts);
```

## Notes

- This framework does not change model physics.
- It reads from `coag_model_final/`.
- It uses the new paper-side PSD diagnostic helpers in `frag_observability/diagnostics/`.
