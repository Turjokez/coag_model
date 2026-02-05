classdef MassBalanceBiovolume
    %MASSBALANCEBIOVOLUME Diagnostics in biovolume-consistent units
    %
    % State vector Y is in "state units" (number per section). This class
    % converts to biovolume using grid.av_vol and computes mass-balance
    % terms in consistent biovolume/time units.

    methods (Static)
        function budget = compute(sim, res, varargin)
            % Optional: state_is_biovolume (default true)
            % If true, Y is already a biovolume/volume state.
            % If false, convert number -> biovolume using grid.av_vol.
            p = inputParser;
            addParameter(p, 'state_is_biovolume', true, @(x)islogical(x) || isnumeric(x));
            parse(p, varargin{:});

            state_is_biovolume = logical(p.Results.state_is_biovolume);

            cfg  = sim.config; %#ok<NASGU>
            grid = sim.grid;

            t = res.time(:);
            Y = res.concentrations;
            n_times = numel(t);
            n_sec   = size(Y, 2);

            if state_is_biovolume
                vbin = ones(n_sec, 1);
            else
                vbin = grid.av_vol(:); % cm^3 per particle (per section)
                if numel(vbin) ~= n_sec
                    error('MassBalanceBiovolume: size mismatch vbin vs Y');
                end
            end

            % Inventory (biovolume-consistent)
            M = Y * vbin; % n_times x 1
            dMdt = gradient(M, t);

            % RHS terms in biovolume units
            rate_coag   = zeros(n_times, 1);
            rate_linear = zeros(n_times, 1);
            rate_disagg = zeros(n_times, 1);
            rate_pp     = zeros(n_times, 1);

            for i = 1:n_times
                [term1, term2, term3, term4, term5] = sim.rhs.rateTerms(Y(i, :)');
                rate_coag(i)   = sum((term1 + term2) .* vbin);
                rate_linear(i) = sum(term3 .* vbin);
                rate_disagg(i) = sum(term4 .* vbin);
                rate_pp(i)     = sum(term5 .* vbin);
            end

            rate_total = rate_coag + rate_linear + rate_disagg + rate_pp;
            residual  = dMdt - rate_total;

            % Integrated (robust) residual: M(t) vs integral of RHS
            M_pred = M(1) + cumtrapz(t, rate_total);
            residual_int = M - M_pred;

            % Midpoint residual (less noisy than gradient)
            dMdt_fd = diff(M) ./ diff(t);
            rate_mid = 0.5 * (rate_total(1:end-1) + rate_total(2:end));
            residual_mid = dMdt_fd - rate_mid;

            % Sinking loss (biovolume/time) using RHS-consistent sinkLossSect
            sink_loss = [];
            if isfield(res, 'output_data') && isfield(res.output_data, 'sinkLossSect')
                sink_loss = res.output_data.sinkLossSect * vbin; % n_times x 1
            elseif isfield(sim, 'operators') && isfield(sim.operators, 'sink_loss')
                sink_diag = diag(sim.operators.sink_loss);
                if isrow(sink_diag), sink_diag = sink_diag'; end
                sink_loss = sum((Y .* sink_diag') .* vbin', 2);
            end

            growth_rate = [];
            if ~isempty(sink_loss)
                % linear = growth - sink  => growth = linear + sink
                growth_rate = rate_linear + sink_loss;
            end

            budget = struct();
            budget.t = t;
            budget.vbin = vbin;
            budget.inventory = M;
            budget.dMdt = dMdt;
            budget.rate_coag = rate_coag;
            budget.rate_linear = rate_linear;
            budget.rate_disagg = rate_disagg;
            budget.rate_pp = rate_pp;
            budget.rate_total = rate_total;
            budget.residual = residual;
            budget.residual_int = residual_int;
            budget.residual_mid = residual_mid;
            budget.rate_mid = rate_mid;
            budget.dMdt_fd = dMdt_fd;
            budget.sink_loss = sink_loss;
            budget.growth_rate = growth_rate;
        end

        function plotSummary(budget, outdir, tag)
            if nargin < 2 || isempty(outdir)
                outdir = fullfile(pwd, 'output', 'figures');
            end
            if nargin < 3
                tag = '';
            end
            if ~isempty(tag)
                tag = ['_' tag];
            end
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end

            t = budget.t;

            % Inventory
            figure('Name', 'Biovolume Inventory');
            plot(t, budget.inventory, 'LineWidth', 1.5);
            grid on;
            xlabel('Time (days)');
            ylabel('Biovolume (state units)');
            title('Total biovolume inventory');
            saveas(gcf, fullfile(outdir, ['mass_balance_inventory' tag '.png']));

            % dM/dt vs total RHS
            figure('Name', 'Biovolume Budget');
            plot(t, budget.dMdt, 'LineWidth', 1.5, 'DisplayName', 'dM/dt');
            hold on;
            plot(t, budget.rate_total, 'LineWidth', 1.5, 'DisplayName', 'sum RHS terms');
            plot(t, budget.residual, 'LineWidth', 1.0, 'DisplayName', 'residual');
            grid on;
            xlabel('Time (days)');
            ylabel('Biovolume rate (state units/day)');
            title('Biovolume mass balance');
            legend('Location', 'best');
            saveas(gcf, fullfile(outdir, ['mass_balance_budget' tag '.png']));

            % Integrated residual (mass closure)
            if isfield(budget, 'residual_int') && ~isempty(budget.residual_int)
                figure('Name', 'Biovolume Budget (Integrated Residual)');
                plot(t, budget.residual_int, 'LineWidth', 1.5);
                grid on;
                xlabel('Time (days)');
                ylabel('Biovolume (state units)');
                title('Integrated mass-balance residual (M - M_{pred})');
                saveas(gcf, fullfile(outdir, ['mass_balance_residual_int' tag '.png']));
            end

            % PP vs losses (sinking + disagg if available)
            if ~isempty(budget.sink_loss)
                show_disagg = isfield(budget, 'rate_disagg') && ~isempty(budget.rate_disagg) ...
                    && any(budget.rate_disagg ~= 0);
                if show_disagg
                    figure('Name', 'PP vs Losses');
                else
                    figure('Name', 'PP vs Sinking Export');
                end
                plot(t, budget.rate_pp, 'LineWidth', 1.5, 'DisplayName', 'PP input');
                hold on;
                plot(t, budget.sink_loss, 'LineWidth', 1.5, 'DisplayName', 'Sinking export');
                if show_disagg
                    disagg_loss = -budget.rate_disagg; % treat disagg as physical loss
                    plot(t, disagg_loss, 'LineWidth', 1.5, 'DisplayName', 'Disagg loss');
                end
                if ~isempty(budget.growth_rate)
                    plot(t, budget.growth_rate, '--', 'LineWidth', 1.0, 'DisplayName', 'Growth (derived)');
                end
                grid on;
                xlabel('Time (days)');
                ylabel('Biovolume rate (state units/day)');
                if show_disagg
                    title('PP input vs losses (sinking + disagg)');
                else
                    title('PP input vs sinking export (biovolume-consistent)');
                end
                legend('Location', 'best');
                saveas(gcf, fullfile(outdir, ['mass_balance_pp_vs_export' tag '.png']));
            end
        end

        function plotRateTermsBySize(sim, res, times, outdir, tag)
            if nargin < 3 || isempty(times)
                times = [0 5 10 20 30 40 50];
            end
            if nargin < 4 || isempty(outdir)
                outdir = fullfile(pwd, 'output', 'figures');
            end
            if nargin < 5
                tag = '';
            end
            if ~isempty(tag)
                tag = ['_' tag];
            end
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end

            t = res.time(:);
            Y = res.concentrations;
            gridObj = sim.grid;
            vbin = gridObj.av_vol(:);

            D_mm = gridObj.getVolumeDiameters() * 10; % cm -> mm

            for k = 1:numel(times)
                [~, idx] = min(abs(t - times(k)));

                [term1, term2, term3, term4, term5] = sim.rhs.rateTerms(Y(idx, :)');
                coag_bin  = abs((term1 + term2) .* vbin);
                lin_bin   = abs(term3 .* vbin);
                dis_bin   = abs(term4 .* vbin);
                pp_bin    = abs(term5 .* vbin);

                % sink loss (if available)
                sink_bin = [];
                if isfield(res, 'output_data') && isfield(res.output_data, 'sinkLossSect')
                    sink_bin = abs(res.output_data.sinkLossSect(idx, :)' .* vbin);
                end

                figure('Name', sprintf('Rate terms by size (t=%g)', t(idx)));
                loglog(D_mm, coag_bin, 'LineWidth', 1.5, 'DisplayName', 'Coag');
                hold on;
                if ~isempty(sink_bin)
                    loglog(D_mm, sink_bin, 'LineWidth', 1.5, 'DisplayName', 'Sink');
                end
                loglog(D_mm, pp_bin, 'LineWidth', 1.5, 'DisplayName', 'PP');
                if any(dis_bin > 0)
                    loglog(D_mm, dis_bin, 'LineWidth', 1.5, 'DisplayName', 'Disagg');
                end
                loglog(D_mm, lin_bin, '--', 'LineWidth', 1.0, 'DisplayName', 'Linear (growth-sink)');

                grid on;
                xlabel('Diameter (mm)');
                ylabel('|dQ/dt| (biovolume units/day)');
                title(sprintf('Rate terms by size (day %g)', t(idx)));
                legend('Location', 'best');

                fname = sprintf('rate_terms_by_size_day_%03d%s.png', round(t(idx)), tag);
                saveas(gcf, fullfile(outdir, fname));
            end
        end
    end
end
