%RUN_ALL_ESTIMATIONS Grid search on the primary cohort, then estimation of
% every analysis defined in config/settings.json.
% Run from MATLAB after R/01_build_sequences.R. Set run_grid = false to skip
% the grid search and use the penalties stored in config/settings.json.

here = fileparts(mfilename('fullpath'));
addpath(here); setup_paths();
cfg = read_settings();

run_grid = true;
if run_grid
    run_grid_search('baseline');
end

analyses = fieldnames(cfg.analyses);
for i = 1:numel(analyses)
    run_estimation(analyses{i});
end
