function folder = campaign_path(model_id,noise_std,generator,category)
%CAMPAIGN_PATH Deterministic directory; no absolute workstation paths.
if nargin<4,category='monte_carlo';end
root=fileparts(fileparts(mfilename('fullpath')));
tag=strrep(sprintf('%.9g',noise_std),'.','p');
folder=fullfile(root,'outputs',category,[model_id '_' generator '_sigma_' tag]);
end
