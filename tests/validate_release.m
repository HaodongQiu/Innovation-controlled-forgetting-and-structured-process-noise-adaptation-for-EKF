function report = validate_release
%VALIDATE_RELEASE Portable syntax/data checks and a short six-method smoke test.
% The shortened test is a diagnostic, not one of the publication's cases.
root=setup_project;files=dir(fullfile(root,'**','*.m'));messages=struct([]);
for k=1:numel(files)
 findings=checkcode(fullfile(files(k).folder,files(k).name),'-id');
 for j=1:numel(findings)
  if contains(findings(j).message,'Parse error') || contains(findings(j).message,'syntax error')
   error('MATLAB:ReleaseSyntax','%s:%d %s',files(k).name,findings(j).line,findings(j).message);
  end
 end
 if ~isempty(findings)
  messages(end+1).file=files(k).name;messages(end).findings=findings; %#ok<AGROW>
 end
end
methods=paper_methods;records=struct([]);
for name={'Pirelli39','ASCE4'}
 model=paper_model(name{1});spec=paper_case(model,1);
 assert(numel(spec.damaged_stories)==(1+2*strcmp(name{1},'Pirelli39')));
 assert(all(spec.true_k>0));
 model.N=600;model.half_N=300;model.forced_samples=450;
 model.ag=model.ag(1:450);model.ag_full=[model.ag zeros(1,150)];
 for k=1:numel(methods)
  cfg=paper_settings(model,.01,methods(k));cfg.windows=[1 600];cfg.window_names={'smoke'};
  [d,v,a,F]=newmark_response(model,spec.true_k);
  [xf,~,nv,h,m,diagnostic]=aekf_filter(model.M,model.K0,d,v,a,F,model.dt,spec.noise_seed,.01,cfg,methods(k),spec.true_k);
  assert(diagnostic.complete && all(isfinite(xf)),'Short test failed: %s %s',name{1},methods(k).id);
  assert(all(nv>0) && m.count==600 && size(h.response_estimate,1)==3*model.n);
  if strcmp(methods(k).controller,'fixed')
   assert(all(h.full_applied_alpha==.8),'Fixed factor differs from .80.');
  else
   assert(all(h.full_applied_alpha>=.8 & h.full_applied_alpha<=.95));
  end
  entry=struct('model',name{1},'method',methods(k).id,'passed',true);
  if isempty(records),records=entry;else,records(end+1)=entry;end %#ok<AGROW>
 end
end
assert(isequal(sample_percentile((1:4)',[0 25 50 75 100]),[1 1.5 2.5 3.5 4]));
report=struct('MATLAB_version',version,'short_tests',records,'code_analyzer_notes',messages, ...
 'full_campaign_rerun',false,'generated',datestr(now,30));
folder=fullfile(root,'outputs','validation');if ~isfolder(folder),mkdir(folder);end
f=fopen(fullfile(folder,'VALIDATION.json'),'w');assert(f>=0);cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('Validation passed: 12 short model/configuration checks.\n');
end
