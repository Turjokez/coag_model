function b = slams_beta(r1, r2, w1, w2, cfg)
% Compute coagulation kernel beta [cm^3/s] for a particle pair.
% Uses same physics as sectional model: Brownian + curvilinear shear + DS.
%
% Inputs:
%   r1, r2  : fractal radii [cm]
%   w1, w2  : settling speeds [cm/s]
%   cfg     : SimulationConfig
%
% Output:
%   b : total kernel [cm^3/s], alpha already applied

    mu = cfg.kvisc * cfg.rho_fl;      % dynamic viscosity [g/(cm s)]

    % Brownian (Smoluchowski)
    conBr = (2/3) * cfg.k * cfg.temp / mu;
    b_br  = conBr * (r1 + r2)^2 / (r1 * r2);

    % Curvilinear shear (Jackson 1990)
    rg  = (r1 + r2) * cfg.r_to_rg;
    p   = min(r1, r2) / max(r1, r2);
    eff = 1 - (1 + 5*p + 2.5*p^2) / (1 + p)^5;
    b_sh = sqrt(8*pi/15) * cfg.gamma * eff * rg^3;

    % Differential settling (curvilinear, small sphere capture)
    r_sm = min(r1, r2) * cfg.r_to_rg;
    b_ds = 0.5 * pi * r_sm^2 * abs(w1 - w2);

    b = cfg.alpha * (b_br + b_sh + b_ds);
end
