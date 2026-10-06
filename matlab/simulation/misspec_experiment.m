function metrics = misspec_experiment(sc, out_dir)
%MISSPEC_EXPERIMENT Sensitivity to the choice of kernel basis (landmark
% positions, number of landmarks and bandwidth).
%   1. calibrate a ground truth with landmarks sc.landmark and bandwidth sc.w
%      (sc.truth = 1: decay 0.5^l as in the main simulation study;
%       sc.truth = 2 or 3: delayed excitation, calibrate_hawkes_params_delayed);
%   2. simulate sc.N sequences; split 60/20/20;
%   3. select (alphaS, alphaGS) by validation BIC with the correctly specified basis;
%   4. for each estimation basis in sc.est (including the matched one), fit on
%      training + validation with the selected penalties and compare with the
%      ground truth. Kernels are compared on the common lag grid 0..M-1, so the
%      comparison is valid when the estimation basis differs from the true one.

if ~isfolder(out_dir), mkdir(out_dir); end
rng(sc.seed);
t0 = tic;

options = struct('N', sc.N, 'Tmax', sc.Tmax, 'dt', 1, 'tstep', 1, 'M', sc.M, 'Nmax', sc.Nmax);
L = numel(sc.landmark);
if sc.truth == 1
    para1 = calibrate_hawkes_params(sc.D, L, 'gauss', sc.w, sc.landmark, options, ...
                                    sc.target_mean, sc.target_sd, sc.attempts, sc.sparsity);
else
    para1 = calibrate_hawkes_params_delayed(sc.D, L, 'gauss', sc.w, sc.landmark, options, ...
                                    sc.target_mean, sc.target_sd, sc.attempts, sc.sparsity, sc.truth);
end
Seqs = Simulation_Thinning_HP_Safe(para1, options);
[A_true, Phi_true] = ImpactFunc(para1, options);

N_train = floor(0.6 * sc.N);
N_val   = floor(0.2 * sc.N);
Seqs_train = Seqs(1:N_train);
Seqs_val   = Seqs(N_train+1:N_train+N_val);
Seqs_test  = Seqs(N_train+N_val+1:end);
Seqs_trval = [Seqs_train, Seqs_val];

alg0 = struct('LowRank', 0, 'Sparse', 1, 'GroupSparse', 1, 'rho', sc.rho, ...
              'outer', sc.outer, 'inner', sc.inner, 'thres', sc.thres, ...
              'Tmax', [], 'storeErr', 0, 'storeLL', 0);

% penalties selected with the correctly specified basis
n_val_events = sum(arrayfun(@(s) numel(s.Time), Seqs_val));
best_bic = Inf;
for aS = sc.alphaS
    for aGS = sc.alphaGS
        alg = alg0; alg.alphaS = aS; alg.alphaGS = aGS;
        m = Initialization_Basis(Seqs_train, 'gauss', sc.w, sc.landmark);
        m = Learning_MLE_Basis(Seqs_train, m, alg);
        bic = -2 * Loglike_Basis(Seqs_val, m, alg) + (sc.D + sc.D^2 * L) * log(n_val_events);
        if bic < best_bic, best_bic = bic; best = [aS, aGS]; end
    end
end
fprintf('Selected penalties: alphaS = %g, alphaGS = %g\n', best(1), best(2));

E_true = A_true > sc.edge_threshold;
rows = cell(numel(sc.est), 1);
for k = 1:numel(sc.est)
    ec = sc.est(k);
    alg = alg0; alg.alphaS = best(1); alg.alphaGS = best(2);
    model = Initialization_Basis(Seqs_trval, 'gauss', ec.w, ec.landmark);
    model = Learning_MLE_Basis(Seqs_trval, model, alg);
    [A_hat, Phi_hat] = ImpactFunc(model, options);

    E_hat = A_hat > sc.edge_threshold;
    TP = sum(E_true(:) & E_hat(:)); FP = sum(~E_true(:) & E_hat(:)); FN = sum(E_true(:) & ~E_hat(:));
    precision = TP / max(TP + FP, 1e-10);
    recall    = TP / max(TP + FN, 1e-10);

    % lag of the kernel peak on the true edges (months)
    [~, pk_true] = max(Phi_true, [], 2);
    [~, pk_hat]  = max(Phi_hat,  [], 2);
    pk_err = abs(squeeze(pk_true) - squeeze(pk_hat));
    pk_err = mean(pk_err(E_true'));                 % Phi is (cause, lag, effect); A is (effect, cause)

    rows{k} = table(string(ec.label), numel(ec.landmark), ec.w, precision, recall, ...
        2 * precision * recall / max(precision + recall, 1e-10), ...
        norm(model.mu - para1.mu) / norm(para1.mu), ...
        norm(Phi_hat(:) - Phi_true(:)) / norm(Phi_true(:)) / sc.D^2, ...
        norm(A_hat(:) - A_true(:)) / norm(A_true(:)), pk_err, ...
        Loglike_Basis(Seqs_test, model, alg), ...
        'VariableNames', {'estimation_basis', 'L', 'w', 'precision', 'recall', 'f1', ...
        'error_mu', 'error_phi', 'error_integrated', 'peak_lag_error_months', 'test_loglik'});
    fprintf('  %-20s F1 = %.3f, error_phi = %.2e\n', ec.label, rows{k}.f1, rows{k}.error_phi);
end

metrics = vertcat(rows{:});
metrics = [table(repmat(string(sc.name), height(metrics), 1), 'VariableNames', {'scenario'}), metrics];
writetable(metrics, fullfile(out_dir, 'metrics.csv'));
save(fullfile(out_dir, 'result.mat'), 'sc', 'para1', 'A_true', 'best', 'metrics');
fprintf('Scenario %s done in %.1f min\n', sc.name, toc(t0) / 60);
end
