%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First add manually:                                                     %
% - MOSEK Solver (https://www.mosek.com/downloads/list/9/)                %
% - SMAC LFR Toolbox v2.1 (https://w3.onera.fr/smac/lfrt)                 %
% - SMAC GSS Library v2.1 (https://w3.onera.fr/smac/gss_download)         %
% to the Matlab path.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath(genpath('/home/ppolcz/_/toolboxes/rsync/01_optimization/mosek'))

ROOT = fileparts(mfilename('fullpath'));
if isempty(ROOT) || contains(ROOT,'/tmp/Editor')
    fname = matlab.desktop.editor.getActiveFilename;
    ROOT = fileparts(fname);
end
cd(ROOT)
setenv('ROOT',ROOT)
addpath(ROOT)

if ~exist('tbxmanager','dir')
    mkdir('tbxmanager')
end

cd tbxmanager
if ~exist('tbxmanager.m','file')
    websave('tbxmanager.m', 'http://www.tbxmanager.com/tbxmanager.m');
    tbxmanager
    savepath
end

if ~exist('toolboxes/yalmip','dir')
    tbxmanager install yalmip
end

if ~exist('toolboxes/sedumi','dir')
    tbxmanager install sedumi
end
cd ..

addpath('tbxmanager')
tbxmanager restorepath

if ~exist('toolboxes','dir')
    mkdir('toolboxes')
end

if ~exist('toolboxes/SMAC_LFR_Toolbox','dir')
    cd toolboxes
    websave SMAC_LFR_Toolbox.zip https://w3.onera.fr/smac/sites/default/files/2023-12/SMAC_LFRT_v2.1.zip
    unzip SMAC_LFR_Toolbox.zip SMAC_LFR_Toolbox
    delete SMAC_LFR_Toolbox.zip
    cd ..
end
if ~exist('toolboxes/SMAC_GSS_Toolbox','dir')
    cd toolboxes
    websave SMAC_GSS_Toolbox https://w3.onera.fr/smac/sites/default/files/2023-12/SMAC_GSS_v2.1.zip
    unzip SMAC_GSS_Toolbox.zip SMAC_GSS_Toolbox
    delete SMAC_GSS_Toolbox.zip
    cd ..
end

addpath(genpath(fullfile(ROOT,'toolboxes')))
addpath([ROOT filesep 'ftools'])
addpath([ROOT filesep 'workspace'])
P_init
G_reset
