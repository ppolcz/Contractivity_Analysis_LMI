%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First add manually:                                                     %
% - MOSEK Solver (https://www.mosek.com/downloads/list/9/)                %
% to the Matlab path.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath(genpath('..../mosek'))

ROOT = fileparts(mfilename('fullpath'));
if isempty(ROOT) || contains(ROOT,'/tmp/Editor')
    fname = matlab.desktop.editor.getActiveFilename;
    ROOT = fileparts(fname);
end
cd(ROOT)
setenv('ROOT',ROOT)
addpath(ROOT)

if ~exist('tbxmanager','dir') || ~exist('tbxmanager/toolboxes','dir') ... and the others

    install

end

addpath('tbxmanager')
tbxmanager restorepath

addpath(genpath(fullfile(ROOT,'toolboxes')))
addpath([ROOT filesep 'ftools'])
addpath([ROOT filesep 'workspace'])

P_init

G_reset
