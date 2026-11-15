function plot_paper_results(folder)
%PLOT_PAPER_RESULTS Generate vector PDFs from saved cases; never runs filters.
setup_project;
if nargin>0
 plot_campaign(folder);return;
end
root=fileparts(mfilename('fullpath'));
campaigns=dir(fullfile(root,'outputs','monte_carlo','*'));
for k=1:numel(campaigns)
 if campaigns(k).isdir && ~ismember(campaigns(k).name,{'.','..'})
  path=fullfile(campaigns(k).folder,campaigns(k).name);
  if ~isempty(result_files(path)),plot_campaign(path);end
 end
end
plot_sensitivity_comparisons;
end
