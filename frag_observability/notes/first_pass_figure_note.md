# First-pass Figure Note

The strongest figure for the paper is currently `fig05_delta_b_delta_kappa_relative.png`.

Why it looks strongest:

- it removes the large baseline offsets among sinking laws
- it shows the fragmentation response as a shift from each law's own no-frag reference
- it makes the image-space contrast easier to read than the raw `(b, kappa)` plot

The strongest single image-space contrast in this run is `kriest_8 | strong` relative to its no-frag baseline.

- baseline image-space `kappa` = `-0.01947`
- case image-space `kappa` = `-0.1013`
- shift in `kappa` = `-0.08186`

Support figures:

- `fig01_b_kappa_separability.png` is still useful as the raw metric-space view.
- `fig02_baseline_psd_fit_residuals.png` is good for showing that slope-only looks too simple.
- `fig03_kappa_vs_fragmentation_strength.png` is good for trend summary.
- `fig04_clean_window_noisy_psd.png` is good for showing what may be lost once the PSD is windowed and sampled.
- `fig05_delta_b_delta_kappa_relative.png` is good for removing law-to-law baseline offsets.
- `fig06_noise_1L_vs_10L_delta_space.png` is good for showing that larger sample volume stabilizes the sampled metric space.

Note:

- this note is written by the paper-figure script
- UVP-like sampling was `true`
- sample volume = `1 L`
- extra sampled delta figure uses `25` repeats at `1` and `10` L
