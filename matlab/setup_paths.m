function setup_paths()
%SETUP_PATHS Add the Hawkes-Process-Toolkit and this repository's MATLAB code
% to the MATLAB path. Run setup/install_toolkit.sh once beforehand.

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
thap = fullfile(repo, 'external', 'Hawkes-Process-Toolkit');

if ~isfolder(thap)
    error(['Hawkes-Process-Toolkit not found in %s.\n' ...
           'Run "bash setup/install_toolkit.sh" from the repository root first.'], thap);
end

% Only the toolkit folders used by this project (estimation, simulation,
% kernels and likelihood); the toolkit's Data/ and Test/ folders are not needed.
for sub = {'BasicFunc', 'Learning', 'Analysis', 'Simulation', 'Visualization'}
    addpath(fullfile(thap, sub{1}));
end
addpath(genpath(here));
end
