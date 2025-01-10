%%
% Created on 2024.10.10. (October 10, Thursday), 15:15

fname = mfilename('fullpath');
if isempty(fname) || contains(fname,'/tmp/Editor')
    fname = matlab.desktop.editor.getActiveFilename;
end
dirname = fileparts(fname);

addpath(fullfile(dirname,'ftools'))

P_init