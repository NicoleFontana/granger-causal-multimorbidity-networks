function res = fit_network(Seqs, est, alphaS, alphaGS)
%FIT_NETWORK Fit a multivariate Hawkes process with Gaussian basis kernels
% by sparse-group-lasso regularised maximum likelihood (Learning_MLE_Basis of
% the Hawkes-Process-Toolkit) and summarise it as a directed network.
%
%   est       estimation settings (cfg.estimation of config/settings.json)
%   alphaS    sparse penalty      (default est.alphaS)
%   alphaGS   group-sparse penalty (default est.alphaGS)
%
% Returns a struct with the fitted model and
%   A    D x D integrated excitation, A(cause, effect)   = edge strengths
%   Phi  D x M x D excitation kernels, Phi(cause, lag, effect), lag in months
%   mu   D x 1 baseline intensities
%   loglik, n_params, n_events, BIC

if nargin < 3, alphaS  = est.alphaS;  end
if nargin < 4, alphaGS = est.alphaGS; end

Tmax     = max([Seqs.Stop]);
landmark = default_landmarks(est.n_landmarks, Tmax);

alg = struct('LowRank', 0, 'Sparse', 1, 'alphaS', alphaS, ...
             'GroupSparse', 1, 'alphaGS', alphaGS, ...
             'outer', est.outer, 'rho', est.rho, 'inner', est.inner, ...
             'thres', est.thres, 'Tmax', [], 'storeErr', 0, 'storeLL', 0);

model = Initialization_Basis(Seqs, est.kernel, est.w, landmark);
model = Learning_MLE_Basis(Seqs, model, alg);

D = size(model.A, 1);
n_events = sum(arrayfun(@(s) numel(s.Time), Seqs));
loglik   = Loglike_Basis(Seqs, model, alg);
n_params = D + D^2 * numel(landmark);

options = struct('Tmax', Tmax, 'M', est.M, 'dt', 1);
[A_effect_cause, Phi] = ImpactFunc(model, options);   % toolkit: A(effect, cause)

res = struct();
res.model    = model;
res.alg      = alg;
res.landmark = landmark;
res.A        = A_effect_cause';                         % A(cause, effect)
res.Phi      = Phi;                                     % Phi(cause, lag, effect)
res.mu       = model.mu(:);
res.loglik   = loglik;
res.n_params = n_params;
res.n_events = n_events;
res.BIC      = -2 * loglik + n_params * log(n_events);
end
