function para1 = calibrate_hawkes_params_delayed(D, L, kernel, w, landmark, options, target_mean, target_std, max_attempts, density_sparsity, variant)
%CALIBRATE_HAWKES_PARAMS_DELAYED Ground truth with delayed excitation, used in
% the landmark-misspecification study. As calibrate_hawkes_params, but:
%   - gentler decay across basis functions (0.8^l instead of 0.5^l);
%   - about 20% of the non-zero pairs get an additional bump at a random late
%     landmark (upper half of the landmarks):
%       variant 2: bump of 30-60% of the first-landmark amplitude;
%       variant 3: bump of 60-90% plus a further 30-60% at the same landmark.

if nargin < 11, variant = 2; end
fprintf('\n=== Calibration (delayed excitation, variant %d): target mean %.0f, SD %.0f ===\n', ...
        variant, target_mean, target_std);

mu_scale = 1.0;
A_scale  = 0.25;
options_test   = options;
options_test.N = 200;
late = ceil(L/2):L;

for attempt = 1:max_attempts
    para1.kernel   = kernel;
    para1.w        = w;
    para1.landmark = landmark;
    para1.mu       = mu_scale * rand(D, 1) / D;

    A_raw = zeros(D, D, L);
    for l = 1:L
        A_raw(:, :, l) = (0.8^l) * (0.5 + ones(D));
    end
    mask  = rand(D) .* double(rand(D) > density_sparsity);
    A_raw = A_raw .* repmat(mask, [1, 1, L]);

    nz  = find(mask > 0);
    sec = nz(randperm(numel(nz), min(round(0.2 * numel(nz)), numel(nz))));
    for p = 1:numel(sec)
        [di, dj] = ind2sub([D, D], sec(p));
        l_late = late(randi(numel(late)));
        primary = A_raw(di, dj, 1);
        if variant == 2
            A_raw(di, dj, l_late) = A_raw(di, dj, l_late) + primary * (0.30 + 0.30 * rand());
        else
            A_raw(di, dj, l_late) = A_raw(di, dj, l_late) + primary * (0.60 + 0.30 * rand());
            A_raw(di, dj, l_late) = A_raw(di, dj, l_late) + primary * (0.30 + 0.30 * rand());
        end
    end

    spec_rad = max(abs(eig(sum(A_raw, 3))));
    if spec_rad < 1e-10
        mu_scale = mu_scale * 1.5;
        continue;
    end
    A_raw = A_scale * A_raw / spec_rad;

    para1.A = zeros(D, L, D);                     % toolkit layout: A(cause, basis, effect)
    for di = 1:D
        for dj = 1:D
            para1.A(dj, :, di) = squeeze(A_raw(di, dj, :));
        end
    end

    Seqs_test = Simulation_Thinning_HP(para1, options_test);
    event_counts = arrayfun(@(s) length(s.Time), Seqs_test);
    if min(event_counts) <= 1
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
