function mb = SectionalMassBalanceDetailed(sim, res)
% SectionalMassBalanceDetailed
% Detailed per-bin in/out terms for diagnostics.
%
% Outputs are [time x section] arrays in model state units/day.

t = res.time(:);
Y = res.concentrations;

n_times = size(Y, 1);
n_sec = size(Y, 2);

betas = sim.operators.betas;
G = sim.operators.growth;
sink_diag = diag(sim.operators.sink_loss);
if isrow(sink_diag), sink_diag = sink_diag'; end

G_pos = max(G, 0);
G_neg = -min(G, 0);

coag_in = zeros(n_times, n_sec);
coag_out = zeros(n_times, n_sec);
growth_in = zeros(n_times, n_sec);
growth_out = zeros(n_times, n_sec);
sink_in = zeros(n_times, n_sec);
sink_out = zeros(n_times, n_sec);
disagg_in = zeros(n_times, n_sec);
disagg_out = zeros(n_times, n_sec);
pp_in = zeros(n_times, n_sec);
pp_out = zeros(n_times, n_sec);
rhs_eval = zeros(n_times, n_sec);

for it = 1:n_times
    v = Y(it, :)';
    v_pos = max(v, eps);
    v_r = v_pos';
    v_shift = [0, v_r(1:n_sec-1)];

    % Coagulation in/out from legacy sectional formulas
    c_in = v_r .* (v_r * betas.b2) + (v_r * betas.b1) .* v_shift;
    c_out = v_r .* (v_r * (betas.b3 + betas.b4 + betas.b5));
    coag_in(it, :) = c_in;
    coag_out(it, :) = c_out;

    % Growth in/out from growth matrix
    growth_in(it, :) = (G_pos * v_pos)';
    growth_out(it, :) = (G_neg * v_pos)';

    % Sinking out from diagonal sink loss
    sink_out(it, :) = (sink_diag .* v_pos)';
    sink_in(it, :) = zeros(1, n_sec);

    % Disagg and PP from RHS terms
    [~, ~, ~, term4, term5] = sim.rhs.rateTerms(v);
    disagg_in(it, :) = max(term4, 0)';
    disagg_out(it, :) = max(-term4, 0)';
    pp_in(it, :) = max(term5, 0)';
    pp_out(it, :) = max(-term5, 0)';

    rhs_eval(it, :) = sim.rhs.evaluate(t(it), v)';
end

total_in = coag_in + growth_in + sink_in + disagg_in + pp_in;
total_out = coag_out + growth_out + sink_out + disagg_out + pp_out;
net_model = total_in - total_out;

% Time derivative from finite difference (for closure check)
dYdt_fd = zeros(n_times, n_sec);
if n_times > 1
    dt = diff(t);
    dYdt_fd(1:end-1, :) = diff(Y, 1, 1) ./ dt;
    dYdt_fd(end, :) = dYdt_fd(end-1, :);
end

closure_fd = dYdt_fd - net_model;
closure_rhs = rhs_eval - net_model;

mb = struct();
mb.t = t;
mb.Y = Y;
mb.coag_in = coag_in;
mb.coag_out = coag_out;
mb.growth_in = growth_in;
mb.growth_out = growth_out;
mb.sink_in = sink_in;
mb.sink_out = sink_out;
mb.disagg_in = disagg_in;
mb.disagg_out = disagg_out;
mb.pp_in = pp_in;
mb.pp_out = pp_out;
mb.total_in = total_in;
mb.total_out = total_out;
mb.net_model = net_model;
mb.dYdt_fd = dYdt_fd;
mb.rhs_eval = rhs_eval;
mb.closure_fd = closure_fd;
mb.closure_rhs = closure_rhs;
end

