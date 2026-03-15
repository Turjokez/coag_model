function ax = plot_b_kappa_separability(fits, ax, show_names)
% plot_b_kappa_separability
% Plot cases in (b, kappa) space.

if nargin < 2 || isempty(ax)
    figure('Color', 'w');
    ax = gca;
end
if nargin < 3
    show_names = false;
end

if isempty(fits)
    return
end

groups = cell(1, numel(fits));
for i = 1:numel(fits)
    groups{i} = get_group(fits(i));
end

[ug, ~, gid] = unique(groups, 'stable');
cols = lines(numel(ug));

hold(ax, 'on');

for g = 1:numel(ug)
    idx = find(gid == g);
    idx = idx(arrayfun(@(k) fits(k).ok, idx));
    if isempty(idx)
        continue
    end

    [~, ord] = sort(arrayfun(@(k) get_level_rank(fits(k)), idx));
    idx = idx(ord);

    col = cols(g,:);
    plot(ax, [fits(idx).b], [fits(idx).kappa], '-', ...
        'Color', col, 'LineWidth', 1.2, 'DisplayName', ug{g});

    for j = 1:numel(idx)
        i = idx(j);
        marker = get_marker(fits(i));
        plot(ax, fits(i).b, fits(i).kappa, marker, ...
            'Color', col, ...
            'MarkerFaceColor', col, ...
            'MarkerSize', 7, ...
            'LineStyle', 'none', ...
            'HandleVisibility', 'off');

        if show_names
            text(ax, fits(i).b, fits(i).kappa, [' ' get_name(fits(i), i)], ...
                'Color', col, 'FontSize', 9);
        end
    end
end

grid(ax, 'on');
xlabel(ax, 'b');
ylabel(ax, 'kappa');
title(ax, '(b, kappa) separability');
legend(ax, 'Location', 'best');
text(ax, 0.02, 0.02, ...
    'marker shape: circle = no-frag, triangle = weak, square = medium, diamond = strong', ...
    'Units', 'normalized', 'FontSize', 8, 'Color', [0.25 0.25 0.25], ...
    'Interpreter', 'none');

end

function name = get_name(fit, idx)
name = sprintf('case %d', idx);
if isfield(fit, 'name') && ~isempty(fit.name)
    name = fit.name;
end
end

function group = get_group(fit)
group = 'all';
if isfield(fit, 'group') && ~isempty(fit.group)
    group = char(string(fit.group));
end
end

function col = get_color(fit, fallback)
col = fallback;
if isfield(fit, 'color') && ~isempty(fit.color)
    col = fit.color;
end
end

function marker = get_marker(fit)
marker = 'o';
if isfield(fit, 'marker') && ~isempty(fit.marker)
    marker = fit.marker;
end
end

function rank = get_level_rank(fit)
rank = 99;
if ~isfield(fit, 'level_name') || isempty(fit.level_name)
    return
end

name = char(string(fit.level_name));
switch lower(name)
    case 'no_frag'
        rank = 1;
    case 'weak'
        rank = 2;
    case 'medium'
        rank = 3;
    case 'strong'
        rank = 4;
end
end
