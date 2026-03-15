function ax = plot_psd_residuals(fits, ax)
% plot_psd_residuals
% Plot residuals and quadratic residual fit.

if nargin < 2 || isempty(ax)
    figure('Color', 'w');
    ax = gca;
end

if isempty(fits)
    return
end

hold(ax, 'on');
cols = lines(numel(fits));

for i = 1:numel(fits)
    if ~fits(i).ok
        continue
    end

    col = get_color(fits(i), cols(i,:));
    name = get_name(fits(i), i);

    plot(ax, fits(i).x, fits(i).resid, '-o', ...
        'Color', col, 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'DisplayName', sprintf('%s | kappa=%.3g', name, fits(i).kappa));
    plot(ax, fits(i).x, fits(i).resid_quad, '--', ...
        'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

yline(ax, 0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid(ax, 'on');
xlabel(ax, 'log10(D)');
ylabel(ax, 'Residual');
title(ax, 'Residuals from fitted power law');
legend(ax, 'Location', 'best');

end

function name = get_name(fit, idx)
name = sprintf('case %d', idx);
if isfield(fit, 'level_name') && ~isempty(fit.level_name)
    name = char(string(fit.level_name));
    return
end
if isfield(fit, 'name') && ~isempty(fit.name)
    name = fit.name;
end
end

function col = get_color(fit, fallback)
col = fallback;
if isfield(fit, 'color') && ~isempty(fit.color)
    col = fit.color;
end
end
