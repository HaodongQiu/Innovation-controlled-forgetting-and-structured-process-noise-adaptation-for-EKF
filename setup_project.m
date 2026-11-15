function root = setup_project
%SETUP_PROJECT Add only the release's source, plotting and test directories.
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'src'),fullfile(root,'plots'),fullfile(root,'tests'));
end
