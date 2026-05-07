classdef ColumnRHS < handle
    % COLUMNRHS  Right-hand side for the 1-D depth column model.
    %
    % One time step:
    %   1. Transport (upwind advection + flux-form diffusion).
    %   2. Process rates (coagulation, disaggregation) at each depth layer.
    %
    % Sinking is transport only — the 0-D coag RHS here has enable_sinking=false.
    %
    % Phase 1: single beta matrix, same process rates at all depths.
    % Phase 2: depth-scaled kernels via brown_scale, shear_scale, ds_scale vectors.
    %
    % Usage:
    %   rhs = ColumnRHS(cfg, size_grid, col_grid, profile);
    %   Y_new = rhs.stepY(Y, dt);

    properties
        cfg         % SimulationConfig (local copy, enable_sinking = false)
        cfg_orig    % original SimulationConfig (with sinking settings intact)
        size_grid   % DerivedGrid
        col_grid    % ColumnGrid
        profile     % DepthProfile
        coag_rhs    % CoagulationRHS — process rates without sinking
        w_z         % n_z x n_sec, sinking speed [m/day]

        % Phase 2: depth-scaling vectors (precomputed, ready to use)
        brown_scale % n_z x 1
        shear_scale % n_z x 1
        ds_scale    % n_z x 1
    end

    methods
        function obj = ColumnRHS(cfg, size_grid, col_grid, profile)
            obj.cfg_orig  = cfg;
            obj.size_grid = size_grid;
            obj.col_grid  = col_grid;
            obj.profile   = profile;

            % local config: disable sinking so transport handles it
            obj.cfg               = cfg.copy();
            obj.cfg.enable_sinking = false;
            obj.cfg.box_depth     = [];

            % build beta matrices (combined: Brownian + shear + DS)
            assembler = BetaAssembler(obj.cfg, size_grid);
            b_brown   = assembler.computeFor('KernelBrown');
            b_shear   = assembler.computeFor('KernelShear');
            b_ds      = assembler.computeFor(cfg.kernel);
            betas     = assembler.combineAndScale(b_brown, b_shear, b_ds);

            % linear matrix (growth only, no sinking)
            lin = LinearProcessBuilder.linearMatrix(obj.cfg, size_grid);

            % disagg matrices (legacy only; operator_split is handled in applyDisaggSplit)
            [Dm, Dp] = LinearProcessBuilder.disaggregationMatrices(obj.cfg);

            obj.coag_rhs = CoagulationRHS(betas, lin, Dm, Dp, obj.cfg, size_grid);

            % sinking speed field: w(k, s) = w_ref(s) * nu_ref / nu(k)
            obj.w_z = obj.buildWindField();

            % depth-scaling vectors (Phase 2 hooks, safe to compute now)
            obj.brown_scale = profile.brownianScale(cfg);
            obj.shear_scale = profile.shearScale(cfg);
            obj.ds_scale    = profile.dsScale(cfg);
        end

        function Y_new = stepY(obj, Y, dt)
            % STEPY  One full time step: transport then process rates.
            % Y: n_z x n_sec
            % dt: day

            % 1. transport (advection + diffusion)
            Kz_z  = obj.profile.Kz;
            Y_new = ColumnTransport.step(Y, obj.w_z, Kz_z, obj.col_grid.dz, dt);

            % 2. process rates at each depth layer
            n_z = obj.col_grid.n_z;
            for k = 1:n_z
                v_k        = Y_new(k, :)';
                dvdt       = obj.coag_rhs.evaluate(0, v_k);
                Y_new(k,:) = max(v_k + dt * dvdt, 0)';
            end

            % 3. operator_split disagg (if enabled)
            if obj.useOperatorSplitDisagg()
                Y_new = obj.applyDisaggSplit(Y_new);
            end
        end

        function w = buildWindField(obj)
            % w(k, s): sinking speed [m/day] at depth k for size bin s.
            % Viscosity correction: w(k) = w_ref * nu_ref / nu(k)
            v_cms  = SettlingVelocityService.velocityForSections(obj.size_grid, obj.cfg_orig);
            w_ref  = (v_cms / 100) * obj.cfg.day_to_sec;  % m/day, n_sec x 1
            nu_ref = obj.cfg_orig.kvisc;                   % cm^2/s
            scale  = nu_ref ./ obj.profile.nu;             % n_z x 1
            w      = scale * w_ref(:)';                    % n_z x n_sec
        end

        function ok = useOperatorSplitDisagg(obj)
            ok = isprop(obj.cfg_orig,'enable_disagg') && obj.cfg_orig.enable_disagg ...
              && isprop(obj.cfg_orig,'disagg_mode') ...
              && strcmpi(obj.cfg_orig.disagg_mode, 'operator_split');
        end

        function Y = applyDisaggSplit(obj, Y)
            % Apply operator_split disagg independently at each depth layer.
            % Phase 1: same D_max at all depths.
            % Phase 2 (TODO): use eps(k) -> D_max(k) per layer.
            n_z = obj.col_grid.n_z;
            for k = 1:n_z
                v_k     = Y(k, :)';
                v_k_new = DisaggregationOperatorSplit.apply(v_k, obj.size_grid, obj.cfg_orig);
                Y(k,:)  = v_k_new';
            end
        end
    end
end
