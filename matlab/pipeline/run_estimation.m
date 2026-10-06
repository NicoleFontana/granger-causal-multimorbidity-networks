function fit = run_estimation(analysis, seq_dir, out_dir)
%RUN_ESTIMATION Fit the network of one analysis with the penalties in
% config/settings.json (estimation.alphaS, estimation.alphaGS; in the paper
% selected by run_grid_search on the primary cohort) on all individuals.
%
% Output in <out_dir>/<analysis>/:
%   model.mat               fitted model, A (cause x effect), Phi, mu, settings
%   excitation_matrix.csv   integrated excitation (edge strengths), rows = cause
%   baseline_intensity.csv  baseline intensity of each event type
%
%   run_estimation('baseline')

cfg = read_settings();
if nargin < 2, seq_dir = fullfile(cfg.repo, 'results', 'sequences'); end
if nargin < 3, out_dir = fullfile(cfg.repo, 'results', 'models');    end

[Seqs, mapping] = load_sequences(fullfile(seq_dir, analysis));
fprintf('Estimating %s: %d individuals, %d event types\n', analysis, numel(Seqs), height(mapping));

tic;
fit = fit_network(Seqs, cfg.estimation);
fit.time_sec   = toc;
fit.analysis   = analysis;
fit.mapping    = mapping;
fit.estimation = cfg.estimation;

names = cellstr(mapping.event_name);
D = size(fit.A, 1);
if numel(names) ~= D
    error('Event mapping has %d types but the model has %d (some types never observed?)', numel(names), D);
end
fprintf('Done in %.1f min: %d of %d directed edges with strength > %g\n', fit.time_sec/60, ...
        nnz(fit.A > cfg.estimation.edge_threshold), D^2, cfg.estimation.edge_threshold);

od = fullfile(out_dir, analysis);
if ~isfolder(od), mkdir(od); end
save(fullfile(od, 'model.mat'), '-struct', 'fit');

A = array2table(fit.A, 'VariableNames', matlab.lang.makeValidName(names), 'RowNames', names);
writetable(A, fullfile(od, 'excitation_matrix.csv'), 'WriteRowNames', true);
writetable(table(mapping.Mark, mapping.event_name, fit.mu, ...
           'VariableNames', {'Mark', 'event_name', 'mu'}), fullfile(od, 'baseline_intensity.csv'));
end
