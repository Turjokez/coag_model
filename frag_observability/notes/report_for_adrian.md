# Short Report For Adrian

## Working question

While I was working on the size-resolved marine particle **box model**, I kept thinking about one question:

**Can PSD shape diagnose fragmentation/disaggregation under sinking-law uncertainty?**

More specifically:

- if fragmentation is active
- but the PSD slope does not change much
- can the PSD **shape** still show it?

To test that, I made a small paper-side analysis inside `frag_observability/` without changing the core model physics.

## What I did

I used the current box model and ran a focused experiment set with:

- 4 fragmentation levels:
  - `no_frag`
  - `weak`
  - `medium`
  - `strong`
- 4 sinking laws:
  - `current`
  - `kriest_8`
  - `kriest_9`
  - `siegel_2025`
- both:
  - image diameter
  - volume diameter
- UVP-like Poisson sampling

I used three simple diagnostics:

- PSD slope `b`
- residuals from the power-law fit
- curvature `kappa` from a quadratic fit to the residuals

The main idea was:

- maybe `b` alone is too weak
- maybe `kappa` carries more process information

## Main result

The first-pass result looks promising.

### 1. Slope `b` alone is not enough

In image space:

- `current` changed almost not at all with fragmentation
- the other 3 laws did change, but not in the same way

So a change in `b` can come from:

- fragmentation
- or the sinking-law choice itself

### 2. `kappa` helps more than `b`

For `kriest_8`, `kriest_9`, and `siegel_2025`, `kappa` changed more clearly than `b`.

Strongest image-space example:

- `kriest_8 no_frag`: `kappa = -0.0195`
- `kriest_8 strong`: `kappa = -0.1013`
- shift: `delta kappa = -0.0819`

This supports the idea that PSD **shape** is more informative than slope alone.

### 3. Volume diameter gives a cleaner signal

The fragmentation signal is stronger in volume diameter than in image diameter.

So the image-diameter conversion seems to weaken the signal, which is important if we want to connect this to observations.

### 4. Noise matters a lot

With UVP-like Poisson sampling:

- `1 L` makes the signal much noisier
- `10 L` is better
- but `10 L` still does not fully recover the clean separation

So the paper should be careful and frame this as an **observability / detectability** result, not a claim that fragmentation is always easy to diagnose.

## My reading of the result

Right now I think the cleanest message is:

**PSD shape helps diagnose fragmentation better than slope alone, but the signal depends on sinking law, diameter definition, and sampling strength.**

I think the strongest figure is the relative-shift plot, because it removes the large baseline differences among sinking laws.

## Best figures

### Main figure

`fig05_delta_b_delta_kappa_relative.png`

Why:

- it removes the large baseline offsets among sinking laws
- it shows the fragmentation response relative to each law's no-frag case
- it is easier to read than the raw `(b, kappa)` plot

![Figure 5](../figures/fig05_delta_b_delta_kappa_relative.png)

### Raw metric-space view

`fig01_b_kappa_separability.png`

Why:

- it shows the full raw `(b, kappa)` space
- it shows that `current` sits apart from the other laws

![Figure 1](../figures/fig01_b_kappa_separability.png)

### Observational realism

`fig06_noise_1L_vs_10L_delta_space.png`

Why:

- it shows that larger sample volume helps
- but also shows that the sampled metric space is still noisy

![Figure 6](../figures/fig06_noise_1L_vs_10L_delta_space.png)

## Where the work stands now

The first-pass analysis is done.

What is already done:

- diagnostics
- controlled experiment set
- summary tables
- first-pass paper figures
- interpretation note

What is still open:

- final paper framing
- figure polishing
- maybe a small observational anchor later

## My suggestion

If this direction makes sense, the next step would be:

1. keep the paper narrow
2. use `b`, residuals, and `kappa` as the main diagnostics
3. use the relative-shift figure as the main result
4. keep the paper framed as a synthetic observability / identifiability study with the box model

## Files

Main results and figures are in:

- `frag_observability/tables/frag_observability_case_summary.csv`
- `frag_observability/figures/fig01_b_kappa_separability.png`
- `frag_observability/figures/fig04_clean_window_noisy_psd.png`
- `frag_observability/figures/fig05_delta_b_delta_kappa_relative.png`
- `frag_observability/figures/fig06_noise_1L_vs_10L_delta_space.png`
