function w = slams_settling(Nn, cfg)
% Settling speed [cm/s] for aggregates using Kriest_8 law.
% Volume-equivalent diameter: d = 2 * a0 * Nn^(1/3)
%
% Inputs:
%   Nn  : number of primary particles per aggregate (scalar or array)
%   cfg : SimulationConfig
%
% Output:
%   w   : settling speed [cm/s]

    a0     = cfg.d0 / 2;
    d_cm   = 2 * a0 * Nn .^ (1/3);       % volume-equivalent diameter [cm]
    w_mday = 66 * d_cm .^ 0.62;           % Kriest_8 [m/day]
    w      = w_mday * 100 / cfg.day_to_sec; % [cm/s]
end
