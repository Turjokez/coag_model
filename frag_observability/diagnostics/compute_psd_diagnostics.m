function fit = compute_psd_diagnostics(D, N, fit_rng)
% compute_psd_diagnostics
% Get PSD slope b, residuals, and quadratic curvature kappa.

fit = struct( ...
    'ok', false, ...
    'fit_range', [], ...
    'n_all', 0, ...
    'n_fit', 0, ...
    'Dall', [], ...
    'Nall', [], ...
    'Dsel', [], ...
    'Nsel', [], ...
    'x', [], ...
    'y', [], ...
    'alpha', NaN, ...
    'b', NaN, ...
    'slope', NaN, ...
    'r2', NaN, ...
    'logN_fit', [], ...
    'Nfit', [], ...
    'Dfit', [], ...
    'resid', [], ...
    'c0', NaN, ...
    'c1', NaN, ...
    'c2', NaN, ...
    'kappa', NaN, ...
    'resid_quad', []);

if nargin < 3 || isempty(fit_rng)
    fit_rng = [-Inf, Inf];
end
fit.fit_range = fit_rng;

D = D(:);
N = N(:);

ok_all = isfinite(D) & isfinite(N) & (D > 0) & (N > 0);
if nnz(ok_all) < 4
    return
end

Dall = D(ok_all);
Nall = N(ok_all);
[Dall, ord_all] = sort(Dall);
Nall = Nall(ord_all);

fit.Dall = Dall;
fit.Nall = Nall;
fit.n_all = numel(Dall);

ok_fit = ok_all & (D >= fit_rng(1)) & (D <= fit_rng(2));
if nnz(ok_fit) < 4
    return
end

Dsel = D(ok_fit);
Nsel = N(ok_fit);
[Dsel, ord_fit] = sort(Dsel);
Nsel = Nsel(ord_fit);

x = log10(Dsel);
y = log10(Nsel);

p = polyfit(x, y, 1);
alpha = p(2);
b = -p(1);
logN_fit = alpha - b * x;
resid = y - logN_fit;

ss_res = sum((y - logN_fit).^2);
ss_tot = sum((y - mean(y)).^2);
if ss_tot > 0
    r2 = 1 - ss_res / ss_tot;
else
    r2 = NaN;
end

q = polyfit(x, resid, 2);
resid_quad = polyval(q, x);

xfit = linspace(min(x), max(x), 100);
Dfit = 10.^xfit;
Nfit = 10.^(alpha - b * xfit);

fit.ok = true;
fit.n_fit = numel(Dsel);
fit.Dsel = Dsel;
fit.Nsel = Nsel;
fit.x = x;
fit.y = y;
fit.alpha = alpha;
fit.b = b;
fit.slope = -b;
fit.r2 = r2;
fit.logN_fit = logN_fit;
fit.resid = resid;
fit.c0 = q(3);
fit.c1 = q(2);
fit.c2 = q(1);
fit.kappa = q(1);
fit.resid_quad = resid_quad;
fit.Dfit = Dfit;
fit.Nfit = Nfit;
