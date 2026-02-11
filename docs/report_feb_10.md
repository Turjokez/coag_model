# Coagulation Model — Feb 10, 2026 

## What I Tested
- Implemented operator‑split disaggregation (ResetInitialCond logic) and kept legacy disaggregation for comparison.
- Compared export fractions and size spectra: legacy vs operator‑split.
- Swept outer timestep (1 hr → 1 day) with a small `Dmax` to force the cap.

## What I Found
- Operator‑split strongly shifts export to small size classes in 0‑D.
- Changing `Dmax` or the outer timestep had only minor impact on export fractions.
- The operator‑split behavior is very different from legacy in 0‑D, so it likely needs evaluation in 1‑D where epsilon varies with depth.

## Figures 

Legacy vs operator‑split export fractions using **image diameter**.  
Operator‑split stays very small‑dominated.  
![Legacy vs operator‑split (image)](figures/export_compare_image_dmax_50.png)

Legacy vs operator‑split export fractions using **volume diameter**.  
Legacy shows large export at late time; operator‑split does not.  
![Legacy vs operator‑split (volume)](figures/export_compare_volume_dmax_50.png)

Mass spectrum at day 100 (volume diameter).  
Operator‑split suppresses the large‑size tail.  
![Mass spectrum comparison](figures/mass_spectrum_compare_dmax_50.png)

Outer‑timestep sweep: **small export fraction (image)**.  
Changing outer timestep does not move the curves much.  
![Outer‑dt sweep small export](figures/disagg_split_outerdt_small_image.png)

Outer‑timestep sweep: **large export fraction (volume)**.  
Large export remains near zero across outer‑dt values.  
![Outer‑dt sweep large export](figures/disagg_split_outerdt_large_volume.png)
