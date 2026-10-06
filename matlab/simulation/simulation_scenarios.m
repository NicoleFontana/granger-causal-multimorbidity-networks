function scenarios = simulation_scenarios(set_name)
%SIMULATION_SCENARIOS Settings of the simulation studies.
%   'smoke'          small scenario to check the code (a few minutes)
%   'paper'          Supplementary Methods S2: D = 29 at imposed sparsity 30%,
%                    50% and 70%, and D = 67 at 30%; N = 500,000 (HPC scale).
%                    The penalties were selected by the full grid search in the
%                    primary scenario (D = 29, 30%) and reused (alphaS =
%                    alphaGS = 100) in the other scenarios.
%   'nd_sweep'       structural recovery as a function of the number of
%                    sequences N and of event types D (17 scenarios, N up to
%                    500,000; penalties fixed at alphaS = alphaGS = 100; one
%                    ground truth per D shared across N)
%   'misspec'        sensitivity to the kernel basis (landmarks, bandwidth):
%                    D = 25, N = 20,000, nine estimation bases, three ground
%                    truths (1: decay 0.5^l; 2 and 3: delayed excitation)
%   'misspec_smoke'  small version of 'misspec'
% Each scenario has a field `type`: 'recovery' (simulation_experiment) or
% 'misspec' (misspec_experiment).

base = struct('name', '', 'type', 'recovery', 'seed', 1, 'D', 29, 'N', 500000, ...
    'Tmax', 480, 'M', 120, 'Nmax', 120, 'target_mean', 10, 'target_sd', 3, 'attempts', 15, ...
    'landmark', [0 6 12 24 36 60 90 120], 'w', 12, 'sparsity', 0.3, ...
    'alphaS', [1 10 50 100], 'alphaGS', [1 10 50 100 250], ...
    'rho', 0.5, 'outer', 8, 'inner', 5, 'thres', 1e-5, 'edge_threshold', 1e-3, ...
    'truth_file', '', 'truth', 1, 'est', []);

switch set_name
    case 'smoke'
        s = base;
        s.name = 'smoke_D5'; s.D = 5; s.N = 2000; s.Nmax = 50; s.target_mean = 8;
        s.alphaS = [10 100]; s.alphaGS = [10 100];
        scenarios = s;

    case 'paper'
        s1 = base; s1.name = 'D29_sparsity30';
        s2 = base; s2.name = 'D29_sparsity50'; s2.sparsity = 0.5; s2.alphaS = 100; s2.alphaGS = 100;
        s3 = base; s3.name = 'D29_sparsity70'; s3.sparsity = 0.7; s3.alphaS = 100; s3.alphaGS = 100;
        s4 = base; s4.name = 'D67_sparsity30'; s4.D = 67; s4.attempts = 100; s4.alphaS = 100; s4.alphaGS = 100;
        scenarios = [s1, s2, s3, s4];

    case 'nd_sweep'
        grid = [10 1000; 10 5000; 10 20000; 10 100000; 10 500000; ...
                25 1000; 25 5000; 25 20000; 25 100000; 25 500000; ...
                50 5000; 50 20000; 50 100000; 50 500000; ...
                75 20000; 75 100000; 75 500000];
        scenarios = repmat(base, 1, size(grid, 1));
        for i = 1:size(grid, 1)
            s = base;
            s.D = grid(i, 1); s.N = grid(i, 2); s.attempts = 20;
            s.alphaS = 100; s.alphaGS = 100;
            s.name = sprintf('ND_D%d_N%d', s.D, s.N);
            s.truth_file = sprintf('truth_D%d.mat', s.D);   % resolved by run_simulations
            scenarios(i) = s;
        end

    case {'misspec', 'misspec_smoke'}
        est = struct('label', {'matched', 'w_narrow', 'w_wide', 'lm_fewer_coarse', ...
                               'lm_fewer_shortterm', 'lm_extra', 'lm_uniform', 'lm_log_short', 'lm_shifted'}, ...
                     'landmark', {[0 6 12 24 36 60 90 120], [0 6 12 24 36 60 90 120], [0 6 12 24 36 60 90 120], ...
                                  [0 12 36 90 120], [0 6 12 24 36 60], ...
                                  [0 3 6 12 18 24 30 36 48 60 75 90 105 120], ...
                                  [0 15 30 45 60 75 90 120], [0 1 2 4 8 16 32 64], [12 24 36 48 60 72 90 120]}, ...
                     'w', {12, 6, 24, 12, 12, 14, 12, 12, 12});
        s = base;
        s.type = 'misspec'; s.D = 25; s.N = 20000; s.Nmax = 80; s.attempts = 20;
        s.alphaS = [10 50 100]; s.alphaGS = [10 50 100]; s.est = est;
        if strcmp(set_name, 'misspec_smoke')
            s.D = 5; s.N = 2000; s.Nmax = 50; s.target_mean = 8;
            s.alphaS = [10 100]; s.alphaGS = 100; s.est = est([1 2 4 8]);
            s.name = 'misspec_smoke_truth2'; s.truth = 2;
            scenarios = s;
        else
            scenarios = repmat(s, 1, 3);
            for t = 1:3
                scenarios(t).truth = t;
                scenarios(t).name = sprintf('misspec_D25_N20000_truth%d', t);
            end
        end

    otherwise
        error('Unknown scenario set "%s".', set_name);
end
end
