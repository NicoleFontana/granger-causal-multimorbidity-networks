function run_simulations(set_name, names)
%RUN_SIMULATIONS Run a set of simulation scenarios (see simulation_scenarios).
%   run_simulations('smoke')                      % quick check
%   run_simulations('paper')                      % Supplementary Methods S2 (HPC scale)
%   run_simulations('paper', {'D29_sparsity30'})  % one scenario
%   run_simulations('nd_sweep')                   % N x D sweep (HPC scale)
%   run_simulations('misspec')                    % kernel-basis misspecification
% Results: results/simulations/<set>/<scenario>/ and
%          results/simulations/<set>/summary.csv
% Completed scenarios (metrics.csv present) are skipped, so a run can be resumed.

if nargin < 1, set_name = 'smoke'; end
here = fileparts(mfilename('fullpath'));
addpath(here); setup_paths();
repo = fileparts(here);

scenarios = simulation_scenarios(set_name);
if nargin >= 2
    scenarios = scenarios(ismember({scenarios.name}, names));
end

out = fullfile(repo, 'results', 'simulations', set_name);
summary = {};
for i = 1:numel(scenarios)
    sc = scenarios(i);
    od = fullfile(out, sc.name);
    if ~isempty(sc.truth_file)
        sc.truth_file = fullfile(out, 'ground_truth', sc.truth_file);
    end
    if isfile(fullfile(od, 'metrics.csv'))
        fprintf('Skipping %s (already completed)\n', sc.name);
        summary{end+1} = readtable(fullfile(od, 'metrics.csv'), 'TextType', 'string'); %#ok<AGROW>
        continue;
    end
    fprintf('\n===== Scenario %s =====\n', sc.name);
    if strcmp(sc.type, 'misspec')
        summary{end+1} = misspec_experiment(sc, od); %#ok<AGROW>
    else
        summary{end+1} = simulation_experiment(sc, od); %#ok<AGROW>
    end
end
writetable(vertcat(summary{:}), fullfile(out, 'summary.csv'));
end
