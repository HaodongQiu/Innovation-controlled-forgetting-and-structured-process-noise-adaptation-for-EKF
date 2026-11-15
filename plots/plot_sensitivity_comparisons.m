function plot_sensitivity_comparisons
%PLOT_SENSITIVITY_COMPARISONS Four noise levels and three truth generators.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'outputs','comparisons');if ~isfolder(folder),mkdir(folder);end
noises=[.005 .01 .03 .05];errors=nan(size(noises));
for k=1:numel(noises)
 path=fullfile(campaign_path('Pirelli39',noises(k),'newmark'),'summary.json');
 if isfile(path)
  s=jsondecode(fileread(path));j=find(strcmp({s.methods.method},'M12-dynamic'));
  if ~isempty(j),errors(k)=s.methods(j).stiffness_nrmse_pct;end
 end
end
if any(isfinite(errors))
 fig=figure('Visible','off','Position',[100 100 640 360]);
 plot(noises,errors,'o-','LineWidth',1.3);xlim([0 .05]);xticks(0:.01:.05);
 xlabel('Measurement-noise standard deviation (m/s^2)');ylabel('Mean final 10 s stiffness NRMSE (%)');grid on
 export_pdf(fig,fullfile(folder,'noise_sensitivity.pdf'));close(fig);
end
rows=struct([]);
for generator={'newmark','ode45','expm_foh'}
 path=fullfile(campaign_path('Pirelli39',.01,generator{1}),'summary.json');
 if ~isfile(path),continue;end
 s=jsondecode(fileread(path));j=find(strcmp({s.methods.method},'M12-dynamic'));
 if isempty(j),continue;end
 row=s.methods(j);row.generator=generator{1};
 if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(folder,'RESPONSE_GENERATOR_COMPARISON.csv'));end
end
