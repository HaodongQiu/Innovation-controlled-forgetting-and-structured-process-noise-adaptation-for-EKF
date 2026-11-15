function run_matched_case(model_id,case_id,noise_std,generator,methods,signature)
%RUN_MATCHED_CASE Generate truth once; all methods receive identical data.
model=paper_model(model_id);spec=paper_case(model,case_id);
folder=campaign_path(model_id,noise_std,generator);
pending=false(1,numel(methods));
for k=1:numel(methods)
 path=fullfile(folder,sprintf('case_%04d_%s.mat',case_id,methods(k).id));
 if isfile(path)
  saved=load(path,'result');r=saved.result;
  assert(strcmp(r.signature,signature) && r.case_id==case_id && ...
   strcmp(r.method.id,methods(k).id) && r.cfg.noise_std==noise_std && ...
   strcmp(r.generator,generator) && strcmp(r.model_id,model_id) && r.completed, ...
   'Saved result belongs to a different experiment. Use a clean output directory.');
 else
  pending(k)=true;
 end
end
if ~any(pending),fprintf('%s case %03d reused.\n',model_id,case_id);return;end
cfg=paper_settings(model,noise_std,methods(1));
[d,v,a,F,info]=paper_truth(model,spec,generator,cfg);
for k=find(pending)
 cfg=paper_settings(model,noise_std,methods(k));
 result=evaluate_filter(model,spec,cfg,methods(k),generator,d,v,a,F,info,signature);
 path=fullfile(folder,sprintf('case_%04d_%s.mat',case_id,methods(k).id));save_result(path,result);
 fprintf('%s case %03d %-12s successful=%d\n',model_id,case_id,methods(k).id,result.success);
end
end
