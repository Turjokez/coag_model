# Diagnostics Summary (0-D)

Date: 2026-02-04
Project: `/Users/turjo/Documents/New project/coag_model_clean`

## Scope
This summary documents 0-D (box) diagnostics and sinking-law sweeps performed to understand why export is dominated by small particles and large sizes rarely appear.

## Baseline configuration used in most runs
- `n_sections = 35`, `t_final = 50 d`, `delta_t = 1 d`
- PP on: `pp_bin = 1`, `pp_source = 1e-8`
- Coag on, sinking on, disagg off, linear on, growth off
- `box_depth = 2000 m`
- Default sinking law: `kriest_8`

## Tests
- All tests pass in `coag_model_clean`.
- PP injection works at `v0 = 0` (fix in `CoagulationRHS.evaluate`).

## Sinking laws tested
- `current`
- `kriest_8`
- `kriest_9`
- `siegel_2025`
- `kriest_8_capped` (cap via `sinking_w_max_mday`)
- `kriest_8_flat` (flat above `sinking_d_flat_cm`)

## Key diagnostics performed (scripts)
- `scripts/run_compare_sinking_laws.m`
- `scripts/run_50day_sinking_sweep.m`
- `scripts/run_export_fraction_timeseries.m`
- `scripts/run_box_depth_sweep.m`
- `scripts/run_sinking_scale_sweep.m`
- `scripts/run_alpha_sweep.m`
- `scripts/run_pp_source_sweep.m`
- `scripts/run_kernel_component_diagnostic.m`
- `scripts/run_kernel_component_sweep.m`
- `scripts/run_r_to_rg_sweep.m`
- `diagnostics/plot_pp_vs_export_biovolume.m`
- `diagnostics/plot_size_distributions_times_image.m`

## Results (high-level)
- **Current law** is far too fast (hundreds of m/day at ~1 mm).
- **Empirical laws (Kriest-8, Kriest-9, Siegel)** reduce speeds to ~15-30 m/day at ~1 mm but **export is still small-dominated**.
- **Timescale ratio** (tau_sink / tau_coag) is **< 1** for ~1-10 mm sizes, even with shear x10, meaning **sinking dominates large bins**.
- **Shear x10** barely changes export fractions.
- **Box depth sweep** and **sinking_scale sweep** show little change in export fractions.
- **PP source sweep** and **alpha sweep** show little change in export fractions.
- **r_to_rg sweep** increases image-based large fraction modestly (still small-dominated).
- **Large initial particles** (placing all mass in bin 10) yields much higher large export fractions, showing the model *can* retain large particles if they exist.
- **Manual sink reduction above 1 mm** yields a large export increase (diagnostic proof that large-size sinking controls export).
- **kriest_8_flat and kriest_8_capped** do not change export because the flattening/cap only affects very large sizes and is too weak.
- **Aggressive flattening at 0.2 mm** still leaves export small-dominated.

## Quantitative highlights
- Baseline `kriest_8` export fractions (image diam): small ~0.97, med ~0.02, large ~0.01
- Shear x10 (image): small ~0.958, med ~0.028, large ~0.014
- Large-initial particles (image): small ~0.482, med ~0.339, large ~0.179
- Manual sink reduction >1 mm x0.1 (image): small ~0.541, med ~0.228, large ~0.231

## Conclusion (0-D)
In a 0-D box model, **large export remains small-dominated even under very slow sinking laws**. The system removes large particles too efficiently. The diagnostics show that **large-size sinking still outpaces coagulation**, so large sizes rarely accumulate.

## Recommended next steps (if/when moving beyond 0-D)
- Re-introduce **column (depth-resolved) sinking** so particles have residence time before loss.
- Use column diagnostics to confirm whether large export increases when loss is distributed by depth.
- Keep 0-D diagnostics for regression tests, but treat 0-D export as a limiting case.

