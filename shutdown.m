%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First add manually:                                                     %
% - MOSEK Solver (https://www.mosek.com/downloads/list/9/)                %
% to the Matlab path.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rmpath(genpath('..../mosek'))

ROOT = fileparts(mfilename('fullpath'));
if isempty(ROOT) || contains(ROOT,'/tmp/Editor')
    fname = matlab.desktop.editor.getActiveFilename;
    ROOT = fileparts(fname);
end

rmpath(genpath(ROOT))
