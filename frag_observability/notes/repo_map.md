# Repo Map

This note maps the files that already control the main parts of the model.

I did not change any core model file.
I only read the repo and made this new `frag_observability/` work area.

## Main split in this repo

- `coag_model_final/` is the main active size-resolved model.
- `adrian_toymodel/` is older layered / EXPORTS-style code with UVP init helpers.

Important:

- `coag_model_final/` is now a size-resolved box model with sinking loss.
- The old layered UVP workflow is mostly in `adrian_toymodel/`.

## 1. Coagulation

These files already control coagulation in the active model:

| File | What it does |
| --- | --- |
| `coag_model_final/src/KernelLibrary.m` | Kernel formulas: Brownian, shear, differential settling. |
| `coag_model_final/src/BetaAssembler.m` | Builds sectional coagulation beta matrices from the kernels. |
| `coag_model_final/src/BetaMatrices.m` | Stores `b1`, `b2`, `b3`, `b4`, `b5`, `b25`. |
| `coag_model_final/src/CoagulationRHS.m` | Uses the beta matrices in `evaluate()` and `rateTerms()`. |
| `coag_model_final/src/CoagulationSimulation.m` | Orchestrates kernel build, RHS setup, and solver run. |

Good starting points:

- `coag_model_final/src/BetaAssembler.m`
- `coag_model_final/src/CoagulationRHS.m`

## 2. Fragmentation / disaggregation

These files already control fragmentation or disaggregation:

| File | What it does |
| --- | --- |
| `coag_model_final/src/SimulationConfig.m` | Main switches and params: `enable_disagg`, `disagg_mode`, `disagg_C`, `disagg_gamma`, `disagg_epsilon`, `disagg_dmax_cm`, `disagg_frac_next`. |
| `coag_model_final/src/Disaggregation.m` | Legacy disaggregation term used inside the RHS. |
| `coag_model_final/src/DisaggregationOperatorSplit.m` | Operator-split redistribution with `Dmax = C epsilon^-gamma` and bin remap. |
| `coag_model_final/src/CoagulationSimulation.m` | Applies the operator-split loop after each outer step. |
| `coag_model_final/src/LinearProcessBuilder.m` | Keeps old matrix-form disagg placeholders `Dminus` and `Dplus`. |

Main diagnostics already present:

- `coag_model_final/diagnostics/run_disagg_operator_split_compare.m`
- `coag_model_final/diagnostics/run_disagg_operator_split_dmax_sweep.m`
- `coag_model_final/diagnostics/run_disagg_operator_split_outerdt_sweep.m`
- `coag_model_final/diagnostics/run_powerlaw_frag_levels.m`

Tests already present:

- `coag_model_final/tests/Test_Disaggregation.m`
- `coag_model_final/tests/Test_DisaggOperatorSplit.m`

## 3. Sinking law

These files already control sinking-law choice and sinking loss:

| File | What it does |
| --- | --- |
| `coag_model_final/src/SimulationConfig.m` | Holds `sinking_law`, `sinking_size`, `sinking_scale`, and other sinking params. |
| `coag_model_final/src/SettlingVelocityService.m` | Main settling-law switch. Supports `current`, `kriest_8`, `kriest_9`, `siegel_2025`, `kriest_8_capped`, `kriest_8_flat`. |
| `coag_model_final/src/LinearProcessBuilder.m` | Converts settling speed into a box-loss sinking matrix. |
| `coag_model_final/src/OutputGenerator.m` | Saves `set_vel`, `sinkLossSect`, `sinkLossTotal`, and flux spectra. |

Main sinking-law diagnostics already present:

- `coag_model_final/diagnostics/run_powerlaw_frag_settling_compare.m`
- `coag_model_final/scripts/run_compare_sinking_laws.m`
- `coag_model_final/scripts/run_50day_sinking_sweep.m`
- `coag_model_final/scripts/run_sinking_scale_sweep.m`

Test already present:

- `coag_model_final/tests/Test_SinkingVelocity.m`

## 4. PSD output

These files already create PSD output and PSD-based diagnostics:

| File | What it does |
| --- | --- |
| `coag_model_final/src/OutputGenerator.m` | Computes `nspec_v`, `nspec_i`, `masspec_v`, `masspec_i`, `fluxspec`, `fluxspec_i`, diameters, and total flux. |
| `coag_model_final/diagnostics/run_powerlaw_frag_levels.m` | Fits PSD slope in log-log space and plots residual shape. |
| `coag_model_final/diagnostics/run_powerlaw_frag_settling_compare.m` | Compares PSD slope across fragmentation and sinking-law cases. |
| `coag_model_final/diagnostics/plot_size_distributions_times.m` | Time slices of PSD in volume diameter. |
| `coag_model_final/diagnostics/plot_size_distributions_times_image.m` | Time slices of PSD in image diameter. |
| `coag_model_final/diagnostics/export_size_classes.m` | Splits flux into small / medium / large size classes. |

Main report files that already describe the PSD result:

- `coag_model_final/docs/report_mar_02.md`
- `coag_model_final/docs/report_mar_03_settling_compare.md`

## 5. Image diameter vs volume diameter

These files already handle the two size definitions:

| File | What it does |
| --- | --- |
| `coag_model_final/src/DerivedGrid.m` | Builds `getImageDiameters()` and `getVolumeDiameters()`. |
| `coag_model_final/src/OutputGenerator.m` | Uses `diam_v`, `diam_i`, and `diaratio` to convert spectra between the two spaces. |
| `coag_model_final/src/SettlingVelocityService.m` | Uses `sinking_size = 'volume'` or `'image'` when applying empirical sinking laws. |
| `coag_model_final/diagnostics/run_powerlaw_frag_levels.m` | Maps `Dmax` from volume-diameter space to nearest image-diameter bin for plotting. |
| `coag_model_final/diagnostics/run_disagg_operator_split_compare.m` | Compares export classes in both image and volume diameter. |

Older UVP conversion helpers in the legacy folder:

- `adrian_toymodel/toymodel/forAB_calcInit/CalcConsConc.m`
- `adrian_toymodel/toymodel/forAB_calcInit/FractalDiamToConsDiam.m`
- `adrian_toymodel/toymodel/forAB_calcInit/CorrectUVPBins.m`

## 6. UVP-like or Poisson-noise sampling

Current synthetic UVP-like / Poisson-style sampling is already in:

| File | What it does |
| --- | --- |
| `coag_model_final/diagnostics/run_uvp_sigma_volume_test.m` | Main UVP-like Poisson test. Uses sample volumes `0.1`, `1`, `10` L and counts bins with `|z| > 2`. |
| `coag_model_final/diagnostics/run_powerlaw_frag_levels.m` | Makes a simple UVP-like Poisson sample and plots sampled residuals. |

Legacy real-UVP input and conversion helpers are in:

| File | What it does |
| --- | --- |
| `adrian_toymodel/toymodel/forAB_calcInit/LoadUVPData.m` | Loads UVP data from `mergedUVPto500m_inEddy_dailyMeans.mat`. |
| `adrian_toymodel/toymodel/forAB_calcInit/CorrectUVPBins.m` | Replaces UVP bin metadata to match model bins. |
| `adrian_toymodel/toymodel/forAB_calcInit/CalcConcentration.m` | Converts `# L^-1 size^-1` to per-bin concentration. |
| `adrian_toymodel/toymodel/forAB_calcInit/CalcConsConc.m` | Maps UVP fractal bins into conserved-volume bins. |
| `adrian_toymodel/toymodel/forAB_calcInit/calcInit.m` | Runs the older UVP init workflow. |

Important:

- I did not find the UVP MAT file itself in this repo.
- So the old UVP path is present in code, but the local dataset is not here now.

## Best entry points for the frag observability paper

If the next step is only new paper-side diagnostics, the cleanest entry files are:

1. `coag_model_final/diagnostics/run_powerlaw_frag_levels.m`
2. `coag_model_final/diagnostics/run_powerlaw_frag_settling_compare.m`
3. `coag_model_final/diagnostics/run_uvp_sigma_volume_test.m`
4. `coag_model_final/src/OutputGenerator.m`
5. `coag_model_final/src/DerivedGrid.m`
6. `coag_model_final/src/DisaggregationOperatorSplit.m`
7. `coag_model_final/src/SettlingVelocityService.m`

## Safe working rule for this new folder

For now:

- keep all new observability work in `frag_observability/`
- read from `coag_model_final/`
- do not edit `coag_model_final/src/` unless needed later
