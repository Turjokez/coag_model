# Initial Findings

This note is based on a fresh run of the paper workflow in `frag_observability/`.

The run used:

- the 4-law core sinking set:
  - `current`
  - `kriest_8`
  - `kriest_9`
  - `siegel_2025`
- 4 fragmentation levels:
  - `no_frag`
  - `weak`
  - `medium`
  - `strong`
- both image diameter and volume diameter
- UVP-like Poisson sampling with `1 L`

This is a focused box-model result summary. It is not a 1-D run and not a full validation study.

## What experiments were run

The full paper-side run did:

- a 3-value sinking-scale scan for each sinking law
- then 4 main cases for each sinking law:
  - `no_frag`
  - `weak`
  - `medium`
  - `strong`
- then diagnostics in:
  - image diameter
  - volume diameter
- then sampled versions with UVP-like Poisson noise

Total model count:

- `4 laws x (3 scale scans + 4 main cases) = 28 model runs`

Main outputs are now saved in:

- `frag_observability/tables/frag_observability_case_summary.csv`
- `frag_observability/tables/frag_observability_baseline_scales.csv`
- `frag_observability/figures/`

## Did slope b change clearly?

Answer: **only for some sinking laws**.

In image diameter:

- `current` stayed almost flat:
  - `4.2625 -> 4.2626`
- `kriest_8` dropped with stronger fragmentation:
  - `3.8218 -> 3.6918`
- `kriest_9` also dropped:
  - `3.8582 -> 3.7417`
- `siegel_2025` also dropped:
  - `3.8307 -> 3.7038`

So `b` does respond in 3 of the 4 laws, but it is not a stable fragmentation indicator across the full 4-law set.

## Did kappa change more clearly?

Answer: **yes, especially outside the `current` law**.

In image diameter:

- `current` stayed almost flat:
  - `-0.7047 -> -0.7048`
- `kriest_8` became much more negative:
  - `-0.0195 -> -0.1013`
- `kriest_9` became more negative:
  - `-0.0299 -> -0.0779`
- `siegel_2025` became more negative:
  - `-0.0240 -> -0.0950`

In volume diameter the signal is stronger:

- `kriest_8`: `-0.0279 -> 0.4367`
- `kriest_9`: `-0.0551 -> 0.4387`
- `siegel_2025`: `-0.0405 -> 0.4374`

So the main result is:

- for the `kriest_8 / kriest_9 / siegel_2025` group, `kappa` changes more clearly than `b`
- for `current`, both `b` and `kappa` are weak in image space until the strongest case

## Did sinking-law uncertainty confuse b?

Answer: **yes**.

The no-frag image-space slopes already differ a lot across laws:

- `current`: `4.2625`
- `kriest_8`: `3.8218`
- `kriest_9`: `3.8582`
- `siegel_2025`: `3.8307`

That means a change in `b` can come from:

- fragmentation strength
- or just switching the sinking law

This is why `b` alone is not enough for the paper question.

## Did (b, kappa) improve separation?

Answer: **yes, but not equally for all laws**.

What looks good:

- in image diameter, the `kriest_8 / kriest_9 / siegel_2025` cases move down in `kappa` as fragmentation gets stronger
- in volume diameter, those same laws separate even more strongly
- `strong` cases are much easier to distinguish than `weak` cases

What still overlaps:

- `weak` and `medium` image-space cases still sit close together
- the `current` law stays far from the others and behaves differently

So the careful conclusion is:

- `(b, kappa)` is better than `b` alone
- but the clean separation is strongest for:
  - stronger fragmentation
  - the `kriest_8 / kriest_9 / siegel_2025` laws
  - volume-diameter space

## Did image-vs-volume or UVP-like noise weaken the signal?

Answer: **yes**.

### image vs volume

Volume diameter gives the cleaner fragmentation signal.

- in image space, `kappa` changes are present but smaller
- in volume space, `kappa` shifts are much larger and easier to separate

So image-diameter conversion does weaken the signal.

### UVP-like / Poisson noise

At `1 L`, the sampled metrics become very noisy.

Examples in sampled image space:

- `current no_frag`: `kappa = 0.7651`
- `current weak`: `kappa = 1.2301`
- `current medium`: `kappa = 1.4241`
- `current strong`: `kappa = 0.2806`

That is much less stable than the clean result.

So the present result is:

- the clean synthetic signal is there
- but `1 L` sampling weakens it a lot
- the paper should probably show:
  - one clean result
  - one observed-window result
  - one noisy result
- and then later test a larger sample volume like `10 L`

## Which figure looks strongest right now?

The strongest paper figure is still:

- `fig05_delta_b_delta_kappa_relative.png`

Why:

- it removes the large law-to-law baseline offsets
- it makes the fragmentation response easier to compare across laws
- it shows more directly that `current` behaves differently from the other 3 laws

The most intuitive support figure is:

- `fig04_clean_window_noisy_psd.png`

Why:

- it shows the clean shape difference
- it shows the residual-window difference
- it shows how much the signal gets weaker after sampling noise

The best new detectability figure is:

- `fig06_noise_1L_vs_10L_delta_space.png`

Why:

- it shows that `10 L` is better than `1 L`
- but it also shows that larger sample volume does not fully recover the clean metric space
- this supports a careful observability framing instead of an overly strong detection claim

## What still needs improvement before this becomes a paper?

Main next steps:

1. Improve the noisy test:
   - keep the new `1 L` vs `10 L` figure
   - but probably add one more larger-volume or multi-repeat summary
   - the `1 L` case is too noisy to carry the paper by itself
2. Make one cleaner separability version:
   - maybe keep the current raw `(b, kappa)` plot
   - and add one second plot with shifts from each law's no-frag baseline
3. Add one simple slope summary figure:
   - `b` vs fragmentation level
   - mostly to show why `b` alone is weaker
4. Decide the main result space:
   - volume space is cleaner
   - image space is closer to observation
   - the paper should likely show both, but not equally
5. Keep the claims careful:
   - this is a focused observability / identifiability study
   - not full validation
   - not a broad BCP paper

## Bottom line

The paper idea still looks good.

The clearest first-pass result is:

- for 3 of the 4 core sinking laws, fragmentation changes `kappa` more clearly than it changes `b`
- sinking-law choice can confuse `b`
- volume-diameter diagnostics are cleaner than image-diameter diagnostics
- `1 L` UVP-like sampling weakens the signal a lot

So the current best paper claim is:

- `PSD shape helps diagnose fragmentation better than slope alone, but the signal depends on sinking law, diameter definition, and sampling strength`
