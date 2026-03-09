# Report - Mar 9, 2026

## what i tested
- Brownian, shear, and differential sedimentation kernel plots.
- compared different settling-speed cases for differential sedimentation.
- grouped the fine spectrum into `N=2`, `N=3`, and `N=4` coarse classes.
- checked conservation for no-frag, one frag case (`eps=1e-6`), and one exact check.

## what is the code to check
- `run_kernel_contour_gallery.m`
- `run_size_class_general_diagnostic.m`
- `size_class_general.m`

## what i found (result)
- Brownian changes more with size than with temperature in this setup.
- Shear gets stronger when partner size and `epsilon` get larger.
- Differential sedimentation changes a lot with the settling-speed case.
- The coarse grouping works well for `N=2`, `N=3`, and `N=4`.
- Amount and coagulation net rate stay conserved at machine-level error (`~1e-21`).
- The exact check gives zero error, and the bad input checks fail correctly.

## figures
Brownian kernel. This shows the kernel changes more with size than with temperature here.
![Brownian kernel](figures/kernel_brownian_temp_size_contours.png)

Shear kernel. This shows the kernel gets stronger when size and `epsilon` get larger.
![Shear kernel](figures/kernel_shear_eps_size_contours.png)

Differential sedimentation comparison. This shows the settling-speed case changes the kernel pattern a lot.
![Differential sedimentation compare](figures/kernel_ds_size_size_settling_laws.png)

Differential sedimentation ratio to the `orig/current` case. This shows where the other settling-speed cases are stronger or weaker.
![Differential sedimentation ratio](figures/kernel_ds_ratio_to_current.png)

Grouped final spectrum for `N=2`, `N=3`, and `N=4`. This shows how the fine bins are merged into larger classes.
![Grouped spectrum](figures/size_class_general_spectrum_groups.png)

Coarse class amount over time. This shows the grouped classes still give a smooth and physical time evolution.
![Coarse amount](figures/size_class_general_amount_timeseries.png)

Coagulation gain, loss, and net rate for a simple grouped case. This shows how material moves into and out of the coarse classes.
![Gain loss net](figures/size_class_general_gain_loss_net.png)

Conservation summary. This shows the amount and rate errors stay at machine level, and the exact check gives zero error.
![Conservation summary](figures/size_class_general_conservation.png)
