function para1 = calibrate_hawkes_params(D, L, kernel, w, landmark, options, target_mean, target_std, max_attempts, density_sparsity)
%CALIBRATE_HAWKES_PARAMS Ground-truth parameters of a multivariate Hawkes
% process whose simulated sequences have about `target_mean` events per
% individual (SD about `target_std`) and fewer than 5% sequences truncated at
% options.Nmax.
%
% Excitation A (D x D x L) decays across basis functions (0.5^l), is masked
% elementwise with probability `density_sparsity` of being zero, has uniform
% random weights and is rescaled so that the spectral radius of sum_l A_l
% equals A_scale (< 1, stationarity). Baseline rates mu and A_scale are
% adjusted iteratively on test simulations of 200 sequences.

if nargin < 9, max_attempts = 10; end

fprintf('\n=== Calibration: target mean %.0f, SD %.0f ===\n', target_mean, target_std);

mu_scale = 1.0;
A_scale  = 0.25;
options_test   = options;
options_test.N = 200;

for attempt = 1:max_attempts
    para1.kernel   = kernel;
    para1.w        = w;
    para1.landmark = landmark;
    para1.mu       = mu_scale * rand(D, 1) / D;

    para1.A = zeros(D, D, L);
    for l = 1:L
        para1.A(:, :, l) = (0.5^l) * (0.5 + ones(D));
    end
    mask = rand(D) .* double(rand(D) > density_sparsity);
    para1.A = para1.A .* repmat(mask, [1, 1, L]);
    para1.A = A_scale * para1.A / max(abs(eig(sum(para1.A, 3))));

    % toolkit layout: A(cause, basis, effect)
    tmp = para1.A;
    para1.A = reshape(para1.A, [D, L, D]);
    for di = 1:D
        for dj = 1:D
            phi = tmp(di, dj, :);
            para1.A(dj, :, di) = phi(:);
        end
    end

    Seqs_test = Simulation_Thinning_HP(para1, options_test);
    event_counts = arrayfun(@(s) length(s.Time), Seqs_test);

    if min(event_counts) <= 1
        fprintf('  Skipped: mu x %.3f too low (sequences with <= 1 event)\n', mu_scale);
        mu_scale = mu_scale * 1.5;
        A_scale  = min(0.5, A_scale * 1.1);
        continue;
    end

    actual_mean   = mean(event_counts);
    actual_std    = std(event_counts);
    saturated_pct = 100 * sum(event_counts == options.Nmax) / length(event_counts);
    fprintf('  Attempt %d: mu x %.3f, A x %.3f -> mean %.1f, SD %.1f, truncated %.0f%%\n', ...
            attempt, mu_scale, A_scale, actual_mean, actual_std, saturated_pct);

    if actual_std  >= target_std * 0.8  && actual_std  <= target_std * 1.2 && ...
       actual_mean >= target_mean * 0.8 && actual_mean <= target_mean * 1.2 && ...
       saturated_pct < 5
        fprintf('  Parameters accepted.\n');
        return;
    end

    if saturated_pct > 20
        mu_scale = mu_scale * 0.6;  A_scale = A_scale * 0.8;
    elseif saturated_pct > 5
        mu_scale = mu_scale * 0.8;  A_scale = A_scale * 0.9;
    elseif actual_mean > target_mean * 1.3
        mu_scale = mu_scale * 0.85;
    elseif actual_mean < target_mean * 0.7
        mu_scale = mu_scale * 1.3;
    elseif actual_std < target_std
        A_scale = A_scale * 0.7;
    end
    mu_scale = max(0.01, min(mu_scale, 5.0));
    A_scale  = max(0.05, min(A_scale, 0.5));
end

warning('Calibration did not meet the targets after %d attempts; using the last parameters.', max_attempts);
end
