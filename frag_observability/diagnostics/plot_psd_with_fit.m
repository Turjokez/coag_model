function ax = plot_psd_with_fit(fits, ax)
% plot_psd_with_fit
% Plot PSD and fitted power-law line.

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

    loglog(ax, fits(i).Dall, fits(i).Nall, '-', ...
        'Color', col, 'LineWidth', 1.4, 'DisplayName', name);
    loglog(ax, fits(i).Dsel, fits(i).Nsel, 'o', ...
        'Color', col, 'MarkerSize', 4, 'HandleVisibility', 'off');
    loglog(ax, fits(i).Dfit, fits(i).Nfit, '--', ...
        'Color', col, 'LineWidth', 1.1, 'HandleVisibility', 'off');
end

grid(ax, 'on');
xlabel(ax, 'Diameter D');
ylabel(ax, 'N(D)');
title(ax, 'PSD with power-law fit');
set(ax, 'XScale', 'log', 'YScale', 'log');
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
