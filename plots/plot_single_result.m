function plot_single_result(result)
%PLOT_SINGLE_RESULT Stiffness ratios, response errors, empirical NIS and alpha.
folder=fullfile(campaign_path(result.model_id,result.cfg.noise_std,result.generator,'single'),'figures');
if ~isfolder(folder),mkdir(folder);end
r=result;t=r.trace.time;
model=paper_model(r.model_id);damaged=r.spec.damaged_stories(:).';
for k=1:numel(damaged)
 story=damaged(k);truth=r.nominal_k(story)*ones(size(t));
 truth(t>=model.first_damaged_sample_time)=r.true_k(story);
 fig=figure('Visible','off','Position',[100 100 760 330]);
 estimate=r.trace.theta(story,:)/r.nominal_k(story);
 sd=r.trace.theta_sd(story,:)/r.nominal_k(story);
 fill([t fliplr(t)],[estimate-2*sd fliplr(estimate+2*sd)], ...
  [.75 .85 .93],'EdgeColor','none');hold on
 plot(t,estimate,'LineWidth',1.1);
 plot(t,truth/r.nominal_k(story),'k--','LineWidth',1.1);
 xlabel('Time (s)');ylabel('Stiffness / healthy stiffness (-)');
 title(sprintf('Story %d',story));legend('Conditional +/-2 SD','Estimated','True','Location','best');
 xlim([0 110]);xticks(0:10:110);grid on
 export_pdf(fig,fullfile(folder,sprintf('case_%04d_story_%02d.pdf',r.case_id,story)));close(fig);
end
fig=figure('Visible','off','Position',[100 100 760 530]);
time=(0:numel(r.trace.full_applied_alpha)-1)*model.dt;
subplot(2,1,1);semilogy(time,r.trace.full_empirical_normalized_nis,'LineWidth',.8);hold on
yline(1,'k--','Expected normalized NIS = 1');
ylim([1e-1 1e5]);yticks(10.^(-1:2:5));
xlim([0 110]);xticks(0:10:110);ylabel('Empirical NIS / n (-)');grid on
subplot(2,1,2);plot(time,r.trace.full_applied_alpha,'LineWidth',1);
xlim([0 110]);xticks(0:10:110);ylim([.8 .95]);yticks(.8:.05:.95);
xlabel('Time (s)');ylabel('Applied forgetting factor (-)');grid on
export_pdf(fig,fullfile(folder,sprintf('case_%04d_NIS_alpha.pdf',r.case_id)));close(fig);
fprintf('Single-case figures: %s\n',folder);
end
