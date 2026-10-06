function cfg = read_settings(repo)
%READ_SETTINGS Read config/settings.json (study design and estimation settings).
if nargin < 1
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
cfg = jsondecode(fileread(fullfile(repo, 'config', 'settings.json')));
cfg.repo = repo;
end
