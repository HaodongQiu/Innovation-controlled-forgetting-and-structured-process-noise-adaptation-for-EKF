function plot_campaign(folder)
%PLOT_CAMPAIGN Comparison CDFs, endpoint ratios and story-wise final NRMSE.
summary=summarize_results(folder);
output=fullfile(folder,'figures');if ~isfolder(output),mkdir(output);end
table=readtable(fullfile(folder,'CASE_METRICS.csv'),'TextType','string');
fig=figure('Visible','off','Position',[100 100 720 390]);hold on
for k=1:numel(summary.methods)
 id=summary.methods(k).method;values=table.stiffness_nrmse_pct(table.method==id & table.success);
 values=sort(values);if isempty(values),continue;end
 stairs([values(1);values(:)],[0;(1:numel(values))'/numel(values)],'LineWidth',1.3,'DisplayName',id);
end
set(gca,'XScale','log');xlabel('Final 10 s stiffness NRMSE (%)');ylabel('Cumulative fraction (-)');
xlim([.01 100]);xticks([.01 .1 1 10 100]);ylim([0 1]);yticks(0:.2:1);
legend('Location','southeast');grid on
export_pdf(fig,fullfile(output,'method_stiffness_CDF.pdf'));close(fig);
for metric={'early_acceleration_nrmse_pct','final_acceleration_nrmse_pct'}
 fig=figure('Visible','off','Position',[100 100 720 390]);hold on
 for k=1:numel(summary.methods)
  id=summary.methods(k).method;
  values=sort(table.(metric{1})(table.method==id & table.success));
  if isempty(values),continue;end
  stairs([values(1);values(:)],[0;(1:numel(values))'/numel(values)],'LineWidth',1.3,'DisplayName',id);
 end
 if startsWith(metric{1},'early'),xlabel('Early acceleration NRMSE (%)');
 else,xlabel('Final 10 s acceleration NRMSE (%)');end
 ylabel('Cumulative fraction (-)');
 ylim([0 1]);yticks(0:.2:1);legend('Location','southeast');grid on
 export_pdf(fig,fullfile(output,[metric{1} '_CDF.pdf']));close(fig);
end
files=result_files(folder);proposed={};
for k=1:numel(files)
 if contains(files(k).name,'_M12-dynamic.mat')
  saved=load(fullfile(folder,files(k).name),'result');
  if saved.result.success,proposed{end+1}=saved.result;end %#ok<AGROW>
 end
end
if isempty(proposed),return;end
n=numel(proposed{1}.nominal_k);endpoint=zeros(n,numel(proposed));errors=endpoint;mask=false(size(endpoint));
for k=1:numel(proposed)
 r=proposed{k};w=find(strcmp(r.cfg.window_names,'final_10_seconds'));
 endpoint(:,k)=r.final_state(3*n+1:end)./r.true_k;
 errors(:,k)=r.metrics.stiffness_window_rmse_pct(:,w);
 mask(r.spec.damaged_stories,k)=true;
end
fig=figure('Visible','off','Position',[100 100 620 400]);
stories=repmat((1:n)',1,numel(proposed));hold on
scatter(stories(~mask),endpoint(~mask),9,[.65 .65 .65],'filled','MarkerFaceAlpha',.25);
scatter(stories(mask),endpoint(mask),16,[0 .4 .65],'filled','MarkerFaceAlpha',.65);yline(1,'k--');
xlabel('Story index');ylabel('Estimated / true stiffness (-)');
xlim([1 n]);xticks(unique([1 5:5:n n]));legend('Unchanged','Damaged','Exact agreement','Location','best');grid on
export_pdf(fig,fullfile(output,'endpoint_stiffness_ratio.pdf'));close(fig);
truth=[];estimates=[];
for k=1:numel(proposed)
 r=proposed{k};i=r.spec.damaged_stories;
 truth=[truth;r.true_k(i)./r.nominal_k(i)]; %#ok<AGROW>
 estimates=[estimates;r.final_state(3*n+i)./r.nominal_k(i)]; %#ok<AGROW>
end
fig=figure('Visible','off','Position',[100 100 500 410]);
scatter(truth,estimates,17,[0 .4 .65],'filled','MarkerFaceAlpha',.65);hold on
plot([.8 1],[.8 1],'k--');xlim([.8 1]);ylim([.8 1]);xticks(.8:.05:1);yticks(.8:.05:1);
xlabel('True / healthy stiffness (-)');ylabel('Estimated / healthy stiffness (-)');grid on
export_pdf(fig,fullfile(output,'damaged_stiffness_ratio.pdf'));close(fig);
fig=figure('Visible','off','Position',[100 100 620 400]);hold on
bands=zeros(n,3);for i=1:n,bands(i,:)=sample_percentile(errors(i,:),[5 50 95]);end
fill([(1:n) fliplr(1:n)],[bands(:,1)' fliplr(bands(:,3)')],[.70 .84 .93],'EdgeColor','none');
plot(1:n,bands(:,2),'Color',[0 .4 .65],'LineWidth',1.4);
xlim([1 n]);xticks(unique([1 5:5:n n]));xlabel('Story index');ylabel('Final 10 s stiffness NRMSE (%)');
legend('5th--95th percentiles','Median','Location','best');grid on
export_pdf(fig,fullfile(output,'story_stiffness_NRMSE.pdf'));close(fig);
fig=figure('Visible','off','Position',[100 100 580 390]);hold on
for group=1:2
 if group==1,values=sort(errors(mask));name='Damaged stories';
 else,values=sort(errors(~mask));name='Unchanged stories';end
 if ~isempty(values),stairs([values(1);values(:)],[0;(1:numel(values))'/numel(values)],'LineWidth',1.3,'DisplayName',name);end
end
set(gca,'XScale','log');xlim([.01 10]);xticks([.01 .1 1 10]);ylim([0 1]);yticks(0:.2:1);
xlabel('Final 10 s stiffness NRMSE (%)');ylabel('Cumulative fraction (-)');legend('Location','southeast');grid on
export_pdf(fig,fullfile(output,'damaged_unchanged_stiffness_CDF.pdf'));close(fig);
for component=1:3
 fig=figure('Visible','off','Position',[100 100 620 390]);hold on
 for window={'early_post_damage','settled_forced','free_decay'}
  values=zeros(1,numel(proposed));
  for k=1:numel(proposed)
   r=proposed{k};w=find(strcmp(r.cfg.window_names,window{1}));
   values(k)=mean(r.metrics.normalized_rmse_pct((component-1)*n+(1:n),w));
  end
  values=sort(values);stairs([values(1);values(:)],[0;(1:numel(values))'/numel(values)],'LineWidth',1.3,'DisplayName',strrep(window{1},'_',' '));
 end
 names={'displacement','velocity','acceleration'};
 xlabel([names{component} ' NRMSE (%)']);ylabel('Cumulative fraction (-)');ylim([0 1]);yticks(0:.2:1);
 legend('Location','southeast');grid on
 export_pdf(fig,fullfile(output,[names{component} '_window_CDF.pdf']));close(fig);
end
plot_violin_errors(proposed,output);
end
