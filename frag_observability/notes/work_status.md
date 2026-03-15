# Work Status

This note is a simple check of where the `frag_observability/` paper work stands now.

## Short answer

For the **first-pass paper plan**, the main work is done.

You now have:

- the paper workspace
- the repo map
- the PSD diagnostics
- the experiment runner
- the summary tables
- the first-pass figure set
- the first findings note
- the report note for Adrian

So you are **past the setup stage**.
You are now in the **interpretation and paper-shaping stage**.

## Done

- `frag_observability/notes/repo_map.md`
- `frag_observability/diagnostics/`
  - `compute_psd_diagnostics.m`
  - `plot_psd_with_fit.m`
  - `plot_psd_residuals.m`
  - `plot_b_kappa_separability.m`
- `frag_observability/experiments/run_frag_observability_framework.m`
- `frag_observability/experiments/generate_first_pass_paper_figures.m`
- tables:
  - `frag_observability_case_summary.csv`
  - `frag_observability_baseline_scales.csv`
  - `frag_observability_case_summary.mat`
- main figures:
  - `fig01_b_kappa_separability.png`
  - `fig02_baseline_psd_fit_residuals.png`
  - `fig03_kappa_vs_fragmentation_strength.png`
  - `fig04_clean_window_noisy_psd.png`
  - `fig05_delta_b_delta_kappa_relative.png`
  - `fig06_noise_1L_vs_10L_delta_space.png`
- notes:
  - `experiment_plan.md`
  - `initial_findings.md`
  - `first_pass_figure_note.md`
  - `how_to_run.md`
  - `report_for_adrian.md`

## Main result status

These parts are now clear enough for a first paper note:

- `b` alone is not enough
- `kappa` is more useful for `kriest_8`, `kriest_9`, and `siegel_2025`
- `current` behaves differently and can hide the signal in image space
- volume diameter gives a cleaner signal than image diameter
- `1 L` noise weakens the signal a lot
- `10 L` is better, but still not perfect

## Still not done

These are not blockers for the first-pass result, but they are still open:

- no observational data anchor yet
- no full paper text yet
- no final paper figure polishing
- no final decision yet on:
  - whether the main text should focus on image space or volume space
  - how much attention `current` should get in the final paper

## What you do not need right now

You do **not** need to:

- change core model physics
- make the model bigger
- turn this into a broad BCP paper
- do a full validation study

## Best next steps

If you want to move this toward a paper, the best next steps are:

1. Show Adrian the short report and figures.
2. Ask if he likes the paper framing:
   - observability / identifiability / detectability
3. If yes, do one small second pass:
   - clean final figure style
   - maybe add one small observed-data anchor later

## Bottom line

You are in a good place.

The first-pass analysis is done.
The next step is not more setup.
The next step is:

- discuss the result
- decide the final framing
- then polish the paper figures and text
