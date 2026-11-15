function plot_violin_errors(results,folder)
%PLOT_VIOLIN_ERRORS Signed errors, 1 s samples from 10 through 103 s.
% All observations enter statistics; both central and full-range PDFs saved.
n=numel(results{1}.nominal_k);pooled=cell(1,4);
for c=1:4,pooled{c}=zeros(n,0);end
for k=1:numel(results)
 r=results{k};time=r.trace.time;
 take=find(time>=10-1e-8 & time<=103+1e-8 & abs(time-round(time))<1e-8);
 assert(numel(take)==94,'Expected 94 integer-second captures per complete case.');
 model=paper_model(r.model_id);
 kt=repmat(r.nominal_k,1,numel(take));kt(:,time(take)>=model.first_damaged_sample_time)=repmat(r.true_k,1,sum(time(take)>=model.first_damaged_sample_time));
 pooled{1}=[pooled{1} 100*(r.trace.theta(:,take)-kt)./r.nominal_k]; %#ok<AGROW>
 error=r.trace.response_estimate(:,take)-r.trace.response_true(:,take);
 for c=1:3,pooled{c+1}=[pooled{c+1} error((c-1)*n+(1:n),:)];end %#ok<AGROW>
end
names={'stiffness','displacement','velocity','acceleration'};
units={'%','m','m/s','m/s^2'};
if n==39
 bandwidths=[.03 7.3e-6 2.1e-5 3.6e-4];
else
 bandwidths=zeros(1,4);
 for c=1:4
  values=pooled{c}(:);spread=min(std(values),diff(sample_percentile(values,[25 75]))/1.349);
  bandwidths(c)=max(.9*spread*numel(values)^(-.2),eps);
 end
end
metadata=struct('sample_times_seconds',10:103,'cases',numel(results), ...
 'observations_per_story',size(pooled{1},2),'bandwidths',bandwidths, ...
 'responses',{names},'units',{units},'kernel','Gaussian','stiffness_denominator','healthy stiffness');
for c=1:4
 X=pooled{c};whiskers=zeros(n,2);stats=zeros(n,3);
 for i=1:n
  x=X(i,:);stats(i,:)=sample_percentile(x,[25 50 75]);iqr=stats(i,3)-stats(i,1);
  whiskers(i,:)=[min(x(x>=stats(i,1)-1.5*iqr)) max(x(x<=stats(i,3)+1.5*iqr))];
 end
 central=max(abs(whiskers),[],'all')*1.08;
 central=max(central,6*bandwidths(c));
 full=max(abs(X),[],'all')*1.02+3*bandwidths(c);
 for view=1:2
  if view==1,span=central;tag='central';else,span=full;tag='full';end
  span=nice_upper(span);grid=linspace(-span,span,501);
  fig=figure('Visible','off','Position',[100 60 700 max(400,24*n)]);hold on
  for i=1:n
   x=X(i,:);density=zeros(size(grid));
   % Batches keep the temporary Gaussian-kernel matrix bounded.
   for start=1:400:numel(x)
    z=(grid(:)-x(start:min(start+399,end)))/bandwidths(c);
    density=density+sum(exp(-.5*z.^2),2)';
   end
   density=density/(numel(x)*bandwidths(c)*sqrt(2*pi));
   width=.39*density/max(density);
   visible=find(density>=max(density)*1e-3);
   visible=visible(1):visible(end);xgrid=grid(visible);shape=width(visible);
   shape([1 end])=0;
   fill([xgrid fliplr(xgrid)],[i+shape i-fliplr(shape)],[.72 .85 .94],'EdgeColor',[.25 .49 .65]);
   plot(whiskers(i,:),[i i],'k-','LineWidth',.8);
   rectangle('Position',[stats(i,1) i-.11 max(stats(i,3)-stats(i,1),eps) .22], ...
    'FaceColor','white','EdgeColor','black','LineWidth',.8);
   plot([stats(i,2) stats(i,2)],[i-.11 i+.11],'k-','LineWidth',1.1);
   plot([whiskers(i,1) whiskers(i,1)],[i-.07 i+.07],'k-');
   plot([whiskers(i,2) whiskers(i,2)],[i-.07 i+.07],'k-');
  end
  xline(0,'k:');xlim([-span span]);xticks(linspace(-span,span,5));
  ylim([.5 n+.5]);yticks(1:n);ylabel('Story index');
  xlabel(sprintf('Signed %s estimation error (%s)',names{c},units{c}));
  export_pdf(fig,fullfile(folder,[names{c} '_KDE_box_' tag '.pdf']));close(fig);
 end
 metadata.([names{c} '_minimum_by_story'])=min(X,[],2);
 metadata.([names{c} '_maximum_by_story'])=max(X,[],2);
end
f=fopen(fullfile(folder,'KDE_METHOD.json'),'w');assert(f>=0);cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(metadata,PrettyPrint=true));
end

function upper = nice_upper(value)
order=10^floor(log10(value));upper=ceil(value/order*2)/2*order;
end
