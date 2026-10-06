function metrics = simulation_experiment(sc, out_dir)
%SIMULATION_EXPERIMENT One simulation scenario (Supplementary Methods S2):
%   1. calibrate ground-truth Hawkes parameters;
%   2. simulate sc.N sequences (Ogata thinning);
%   3. split into training (60%), validation (20%) and test (20%) sets;
%   4. select (alphaS, alphaGS) by the validation BIC over sc.alphaS x sc.alphaGS;
%   5. refit on training + validation with the selected penalties;
%   6. compare the estimate with the ground truth.
%
% Metrics (written to <out_dir>/metrics.csv):
%   precision, recall, f1  edge recovery; an edge is present when its
%                          integrated excitation exceeds sc.edge_threshold
%   true/learned_density   % of the D^2 directed pairs that are edges
%   true/learned_sparsity  100 - density (as reported in the paper)
%   error_mu               ||mu_hat - mu|| / ||mu||
%   error_phi              ||Phi_hat - Phi||_F / ||Phi||_F / D^2 (kernel values
%                          on the lag grid, as in the original code)

if ~isfolder(out_dir), mkdir(out_dir); end
rng(sc.seed);
t0 = tic;

options = struct('N', sc.N, 'Tmax', sc.Tmax, 'dt', 1, 'tstep', 1, 'M', sc.M, 'Nmax', sc.Nmax);
L = numel(sc.landmark);

% 1-2. ground truth and sequences
if isfield(sc, 'truth_file') && ~isempty(sc.truth_file) && isfile(sc.truth_file)
    load(sc.truth_file, 'para1');                 % ground truth shared across scenarios (N x D sweep)
else
    para1 = calibrate_hawkes_params(sc.D, L, 'gauss', sc.w, sc.landmark, options, ...
                                    sc.target_mean, sc.target_sd, sc.attempts, sc.sparsity);
    if isfield(sc, 'truth_file') && ~isempty(sc.truth_file)
        if ~isfolder(fileparts(sc.truth_file)), mkdir(fileparts(sc.truth_file)); end
        save(sc.truth_file, 'para1');
    end
end
Seqs = Simulation_Thinning_HP_Safe(para1, options);
[A_true, Phi_true] = ImpactFunc(para1, options);

% 3. split
N_train = floor(0.6 * sc.N);
N_val   = floor(0.2 * sc.N);
Seqs_train = Seqs(1:N_train);
Seqs_val   = Seqs(N_train+1:N_train+N_val);
Seqs_test  = Seqs(N_train+N_val+1:end);

alg0 = struct('LowRank', 0, 'Sparse', 1, 'GroupSparse', 1, 'rho', sc.rho, ...
              'outer', sc.outer, 'inner', sc.inner, 'thres', sc.thres, ...
              'Tmax', [], 'storeErr', 0, 'storeLL', 0);
n_params = sc.D + sc.D^2 * L;
n_val_events = sum(arrayfun(@(s) numel(s.Time), Seqs_val));

% 4. grid search on the validation set
[aS, aGS] = ndgrid(sc.alphaS, sc.alphaGS);
grid = table(aS(:), aGS(:), nan(numel(aS), 1), 'VariableNames', {'alphaS', 'alphaGS', 'val_BIC'});
for i = 1:height(grid)
    alg = alg0; alg.alphaS = grid.alphaS(i); alg.alphaGS = grid.alphaGS(i);
    model = Initialization_Basis(Seqs_train, 'gauss', sc.w, sc.landmark);
    model = Learning_MLE_Basis(Seqs_train, model, alg);
    grid.val_BIC(i) = -2 * Loglike_Basis(Seqs_val, model, alg) + n_params * log(n_val_events);
    fprintf('  alphaS = %g, alphaGS = %g: validation BIC = %.1f\n', grid.alphaS(i), grid.alphaGS(i), grid.val_BIC(i));
end
[~, best] = min(grid.val_BIC);
writetable(grid, fullfile(out_dir, 'grid_search.csv'));

% 5. final fit on training + validation
alg = alg0; alg.alphaS = grid.alphaS(best); alg.alphaGS = grid.alphaGS(best);
Seqs_trval = [Seqs_train, Seqs_val];
model = Initialization_Basis(Seqs_trval, 'gauss', sc.w, sc.landmark);
model = Learning_MLE_Basis(Seqs_trval, model, alg);
[A_hat, Phi_hat] = ImpactFunc(model, options);

% 6. comparison with the ground truth
E_true = A_true > sc.edge_threshold;
E_hat  = A_hat  > sc.edge_threshold;
TP = sum(E_true(:) & E_hat(:)); FP = sum(~E_true(:) & E_hat(:)); FN = sum(E_true(:) & ~E_hat(:));
precision = TP / max(TP + FP, 1e-10);
recall    = TP / max(TP + FN, 1e-10);
n_events  = arrayfun(@(s) numel(s.Time), Seqs);

metrics = table(string(sc.name), sc.D, sc.N, sc.sparsity, ...
    mean(n_events), std(n_events), grid.alphaS(best), grid.alphaGS(best), height(grid), ...
    precision, recall, 2 * precision * recall / max(precision + recall, 1e-10), ...
    100 * mean(E_true(:)), 100 * mean(E_hat(:)), 100 - 100 * mean(E_true(:)), 100 - 100 * mean(E_hat(:)), ...
    norm(model.mu - para1.mu) / norm(para1.mu), ...
    norm(Phi_hat(:) - Phi_true(:)) / norm(Phi_true(:)) / sc.D^2, ...
    Loglike_Basis(Seqs_test, model, alg), toc(t0), ...
    'VariableNames', {'scenario', 'D', 'N', 'imposed_sparsity', 'events_mean', 'events_sd', ...
    'best_alphaS', 'best_alphaGS', 'n_grid_combinations', 'precision', 'recall', 'f1', ...
    'true_density', 'learned_density', 'true_sparsity', 'learned_sparsity', ...
    'error_mu', 'error_phi', 'test_loglik', 'time_sec'});
writetable(metrics, fullfile(out_dir, 'metrics.csv'));
save(fullfile(out_dir, 'result.mat'), 'sc', 'para1', 'model', 'A_true', 'A_hat', 'grid', 'metrics');
disp(metrics);
end
