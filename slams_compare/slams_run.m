function [t_days, Nn_arr, n_mat] = slams_run(cfg, grid, n0, opts)
% SLAMS-style super-particle coagulation model (0-D slab, no sinking).
%
% One super-particle per section. Each particle tracks Nn (primary particles
% per aggregate) and n [#/cm^3] (real aggregate concentration).
% Default fractal dimension is D=2 (SLAMS hardcoded value).
%
% When two types collide, merged particles snap to the nearest existing Nn bin.
% This is the key algorithmic difference from the sectional ODE:
% collisions are counted pair-wise and redistributed discretely.
%
% Inputs:
%   cfg  : SimulationConfig
%   grid : DerivedGrid
%   n0   : initial concentrations [n_sec x 1, #/cm^3]
%   opts : optional struct
%          opts.D         fractal dimension for collision geometry (default 2.0)
%          opts.n_sub     sub-steps per output interval (default 24)
%          opts.snap_mode 'nearest' or 'split' (default 'nearest')
%
% Outputs:
%   t_days : time vector [days]
%   Nn_arr : Nn value for each super-particle slot [n_sec x 1]
%   n_mat  : concentrations [n_times x n_sec, #/cm^3]

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'D') || isempty(opts.D)
        opts.D = 2.0;
    end
    if ~isfield(opts, 'n_sub') || isempty(opts.n_sub)
        opts.n_sub = 24;
    end
    if ~isfield(opts, 'snap_mode') || isempty(opts.snap_mode)
        opts.snap_mode = 'nearest';
    end

    n_sec  = cfg.n_sections;
    a0     = grid.a0;
    D      = opts.D;

    % Super-particle sizes: Nn(k) = 2^(k-1), same doubling as sectional bins
    Nn    = 2 .^ (0:n_sec-1)';             % [n_sec x 1]
    r_agg = a0 * Nn .^ (1/D);              % fractal radius [cm]
    w_agg = slams_settling(Nn, cfg);        % settling speed [cm/s]

    % Pre-compute beta for all pairs [cm^3/s]
    B = zeros(n_sec, n_sec);
    for i = 1:n_sec
        for j = i:n_sec
            bij = slams_beta(r_agg(i), r_agg(j), w_agg(i), w_agg(j), cfg);
            B(i,j) = bij;
            B(j,i) = bij;
        end
    end

    % For each pair (i,j), find the target bin for Nn(i)+Nn(j)
    log_Nn = log(Nn);
    target  = zeros(n_sec, n_sec);
    tlo     = zeros(n_sec, n_sec);
    thi     = zeros(n_sec, n_sec);
    wlo     = zeros(n_sec, n_sec);
    whi     = zeros(n_sec, n_sec);
    for i = 1:n_sec
        for j = i:n_sec
            Nn_new      = Nn(i) + Nn(j);
            [~, k] = min(abs(log_Nn - log(Nn_new)));
            target(i,j) = k; target(j,i) = k;

            if strcmpi(opts.snap_mode, 'split')
                lnew = log(Nn_new);
                if lnew <= log_Nn(1)
                    lo = 1; hi = 1; wl = 1; wh = 0;
                elseif lnew >= log_Nn(end)
                    lo = n_sec; hi = n_sec; wl = 1; wh = 0;
                else
                    hi = find(log_Nn >= lnew, 1, 'first');
                    lo = hi - 1;
                    span = log_Nn(hi) - log_Nn(lo);
                    wh = (lnew - log_Nn(lo)) / span;
                    wl = 1 - wh;
                end
                tlo(i,j) = lo; tlo(j,i) = lo;
                thi(i,j) = hi; thi(j,i) = hi;
                wlo(i,j) = wl; wlo(j,i) = wl;
                whi(i,j) = wh; whi(j,i) = wh;
            end
        end
    end

    % Time grid and sub-stepping
    t_days  = (cfg.t_init : cfg.delta_t : cfg.t_final)';
    n_times = length(t_days);
    n_mat   = zeros(n_times, n_sec);

    n_sub  = opts.n_sub;                    % sub-steps per output interval
    dt_sub = cfg.delta_t / n_sub;           % [days]
    B_day  = B * cfg.day_to_sec;            % [cm^3/day]

    n_now        = n0(:);
    n_mat(1, :)  = n_now';

    for it = 2:n_times
        for sub = 1:n_sub
            dn = zeros(n_sec, 1);

            for i = 1:n_sec
                for j = i:n_sec
                    if n_now(i) <= 0 || n_now(j) <= 0
                        continue
                    end

                    % Merger rate [#/cm^3/day]
                    if i == j
                        Q = 0.5 * B_day(i,j) * n_now(i)^2;
                    else
                        Q = B_day(i,j) * n_now(i) * n_now(j);
                    end

                    dQ = Q * dt_sub;  % mergers this sub-step [#/cm^3]

                    % Limit so concentrations stay non-negative
                    if i == j
                        dQ = min(dQ, 0.5 * n_now(i));
                    else
                        dQ = min(dQ, min(n_now(i), n_now(j)));
                    end

                    % Losses: each merger removes one aggregate from each parent
                    % (for i==j this decrements i twice, which is correct)
                    dn(i) = dn(i) - dQ;
                    dn(j) = dn(j) - dQ;

                    % Gain with mass conservation under snap mapping.
                    Nn_new = Nn(i) + Nn(j);
                    if strcmpi(opts.snap_mode, 'split')
                        lo = tlo(i,j); hi = thi(i,j);
                        wl = wlo(i,j); wh = whi(i,j);
                        dQ_lo = dQ * wl * (Nn_new / Nn(lo));
                        dQ_hi = dQ * wh * (Nn_new / Nn(hi));
                        dn(lo) = dn(lo) + dQ_lo;
                        dn(hi) = dn(hi) + dQ_hi;
                    else
                        k          = target(i,j);
                        mass_scale = Nn_new / Nn(k);
                        dn(k)      = dn(k) + dQ * mass_scale;
                    end
                end
            end

            n_now = max(n_now + dn, 0);
        end

        n_mat(it, :) = n_now';
    end

    Nn_arr = Nn;
end
