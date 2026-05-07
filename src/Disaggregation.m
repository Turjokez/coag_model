classdef Disaggregation
    %DISAGGREGATION Legacy turbulence-driven breakup term (non-mass-conserving).
    %
    % IMPORTANT: This term does NOT conserve mass/biovolume.
    %
    % The formula per bin is:
    %   term(isec) = -c3 * c4^isec * (v(isec) - c4*v(isec+1))   isec = 2..n-1
    %
    % When v(isec) > c4*v(isec+1): term is negative (loss from bin isec).
    % When v(isec) < c4*v(isec+1): term is POSITIVE (creation in bin isec).
    % In neither case is the removed/added mass redistributed to another bin.
    % Mass is not conserved — the term acts as a pure source or sink per bin.
    %
    % Consequence: with c4=1.45 the rate c3*c4^isec grows rapidly with bin
    % index (e.g. c3=0.02, c4=1.45: rate at bin 15 is ~2 day^{-1} while
    % sinking rate at that bin is ~0.015 day^{-1}). This imbalance causes
    % large total-number growth (+800%) at default parameters.
    %
    % Recommended alternatives:
    %   - Reduce c3 until frag rate is comparable to sinking rate per bin.
    %     Use run_may06_frag_rate_analysis.m to diagnose the balance.
    %   - Switch to operator_split disagg (DisaggregationOperatorSplit) which
    %     is mass-conserving and physically based on a D_max threshold.
    %
    % Inputs:
    %   v   - state vector (n x 1)
    %   cfg - SimulationConfig (uses cfg.c3, cfg.c4)

    methods (Static)
        function term = netTerm(v, cfg)
            n = length(v);
            term = zeros(n,1);

            if n < 3
                return;
            end

            c3 = 0.02;
            c4 = 1.45;
            if nargin >= 2 && ~isempty(cfg)
                if isprop(cfg,'c3') && ~isempty(cfg.c3), c3 = cfg.c3; end
                if isprop(cfg,'c4') && ~isempty(cfg.c4), c4 = cfg.c4; end
            end

            % Clip negatives to zero for physics (no EPS injection)
            v_pos = max(v, 0);

            for isec = 2:(n-1)
                term(isec) = term(isec) - c3 * c4^isec * (v_pos(isec) - c4 * v_pos(isec+1));
            end

            % Runtime check: warn if legacy disagg is creating mass this step.
            % net positive sum means the term is adding biovolume (non-physical source).
            net_mass = sum(term);
            if net_mass > 0
                persistent warned_creation;
                if isempty(warned_creation)
                    warned_creation = false;
                end
                if ~warned_creation
                    warning('Disaggregation:massCreation', ...
                        ['Legacy disagg is net-positive this step (net=%.3e). ' ...
                         'Mass is being created, not just redistributed. ' ...
                         'Consider reducing c3 or switching to operator_split mode.'], net_mass);
                    warned_creation = true;
                end
            end
        end
    end
end
