function result = run_single_filter(model_id,method_id,case_id,noise_std,generator)
%RUN_SINGLE_FILTER One complete benchmark case, without generating plots.
% result=run_single_filter('Pirelli39','M12-dynamic',51,0.01,'newmark')
setup_project;
if nargin<1,model_id='Pirelli39';end
if nargin<2,method_id='M12-dynamic';end
if nargin<3,case_id=51;end
if nargin<4,noise_std=.01;end
if nargin<5,generator='newmark';end
model=paper_model(model_id);methods=paper_methods;
take=find(strcmp({methods.id},method_id));assert(isscalar(take),'Unknown method ID. See paper_methods.');
method=methods(take);spec=paper_case(model,case_id);cfg=paper_settings(model,noise_std,method);
[d,v,a,F,info]=paper_truth(model,spec,generator,cfg);
result=evaluate_filter(model,spec,cfg,method,generator,d,v,a,F,info,experiment_signature);
folder=campaign_path(model_id,noise_std,generator,'single');
path=fullfile(folder,sprintf('case_%04d_%s.mat',case_id,method_id));save_result(path,result);
fprintf('%s / %s / case %d: successful=%d. Saved %s\n',model_id,method_id,case_id,result.success,path);
end
