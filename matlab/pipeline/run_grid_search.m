function results = run_grid_search(analysis, seq_dir, out_dir)
%RUN_GRID_SEARCH BIC-based selection of the penalties (alphaS, alphaGS).
% Individuals are split into training and validation sets; each sampled
% combination is fitted on the training set and scored by the BIC computed on
% the validation set, BIC = -2 loglik_val + k log(n_val_events), with
% k = D + D^2 L. Grid, split and sampling are set in config/settings.json.
%
%   run_grid_search('baseline')                       % default folders
%   run_grid_search('baseline', 'results/sequences', 'results/grid_search')

cfg = read_settings();
if nargin < 2, seq_dir = fullfile(cfg.repo, 'results', 'sequences');   end
if nargin < 3, out_dir = fullfile(cfg.repo, 'results', 'grid_search'); end
gs = cfg.grid_search;

Seqs = load_sequences(fullfile(seq_dir, analysis));

rng(gs.seed);
N       = numel(Seqs);
idx     = randperm(N);
n_train = round(gs.train_fraction * N);
Seqs_train = Seqs(idx(1:n_train));
Seqs_val   = Seqs(idx(n_train+1:end));

[aS, aGS] = ndgrid(gs.alphaS, gs.alphaGS);
grid = [aS(:), aGS(:)];
n_comb = size(grid, 1);
n_eval = round(gs.sampling_fraction * n_comb);
rng(gs.seed);
grid = grid(sort(randperm(n_comb, n_eval)), :);

n_val_events = sum(arrayfun(@(s) numel(s.Time), Seqs_val));
results = table(grid(:,1), grid(:,2), nan(n_eval,1), nan(n_eval,1), nan(n_eval,1), ...
    'VariableNames', {'alphaS', 'alphaGS', 'train_loglik', 'val_loglik', 'val_BIC'});

fprintf('Grid search (%s): %d of %d combinations, %d training / %d validation individuals\n', ...
        analysis, n_eval, n_comb, numel(Seqs_train), numel(Seqs_val));
for i = 1:n_eval
    fit = fit_network(Seqs_train, cfg.estimation, grid(i,1), grid(i,2));
    results.train_loglik(i) = fit.loglik;
    results.val_loglik(i)   = Loglike_Basis(Seqs_val, fit.model, fit.alg);
    results.val_BIC(i)      = -2 * results.val_loglik(i) + fit.n_params * log(n_val_events);
    fprintf('  alphaS = %g, alphaGS = %g: validation BIC = %.1f\n', grid(i,1), grid(i,2), results.val_BIC(i));
end

[~, best] = min(results.val_BIC);
fprintf('Selected: alphaS = %g, alphaGS = %g\n', results.alphaS(best), results.alphaGS(best));

if ~isfolder(fullfile(out_dir, analysis)), mkdir(fullfile(out_dir, analysis)); end
writetable(results, fullfile(out_dir, analysis, 'grid_search_results.csv'));
end
