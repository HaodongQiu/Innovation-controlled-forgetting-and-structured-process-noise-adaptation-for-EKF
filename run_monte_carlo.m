function run_monte_carlo(workers,models,noise_std,generator,case_ids,method_ids)
%RUN_MONTE_CARLO Matched 100-case analyses. Resume with the identical command.
% run_monte_carlo(4) runs both structures and all six paper configurations.
setup_project;
if nargin<1,workers=4;end
if nargin<2,models={'Pirelli39','ASCE4'};end
if nargin<3,noise_std=.01;end
if nargin<4,generator='newmark';end
if nargin<5,case_ids=1:100;end
methods=paper_methods;
if nargin<6,method_ids={methods.id};end
if ischar(models)||isstring(models),models=cellstr(models);end
if ischar(method_ids)||isstring(method_ids),method_ids=cellstr(method_ids);end
validateattributes(workers,{'numeric'},{'scalar','integer','positive'});
validateattributes(case_ids,{'numeric'},{'vector','integer','>=',1,'<=',100});
assert(numel(unique(case_ids))==numel(case_ids),'Duplicate case IDs.');
assert(all(ismember(method_ids,{methods.id})),'Unknown method ID.');
methods=methods(ismember({methods.id},method_ids));assert(~isempty(methods));
assert(ismember(generator,{'newmark','ode45','expm_foh'}),'Unknown truth generator.');
if workers>1
 products=ver;
 assert(any(strcmp({products.Name},'Parallel Computing Toolbox')) && ...
  license('test','Distrib_Computing_Toolbox'), ...
  'Parallel Computing Toolbox must be installed and licensed. Use workers=1 for serial execution.');
 pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',workers);end
 assert(pool.NumWorkers==workers,'Existing pool has a different size. Use its count or close it first.');
end
signature=experiment_signature;
for b=1:numel(models)
 model_id=models{b};paper_model(model_id);folder=campaign_path(model_id,noise_std,generator);
 if ~isfolder(folder),mkdir(folder);end
 status=struct('state','running','model_id',model_id,'case_ids',case_ids, ...
  'method_ids',{ {methods.id} },'expected',numel(case_ids)*numel(methods), ...
  'workers',workers,'MATLAB_version',version,'signature',signature);
 write_status(folder,status);
 try
  if workers>1
   parfor k=1:numel(case_ids)
    run_matched_case(model_id,case_ids(k),noise_std,generator,methods,signature);
   end
  else
   for k=1:numel(case_ids)
    run_matched_case(model_id,case_ids(k),noise_std,generator,methods,signature);
   end
  end
  summarize_results(folder);
  status.state='complete';write_status(folder,status);
 catch ME
  status.state='needs_attention';status.error=getReport(ME,'extended','hyperlinks','off');
  write_status(folder,status);rethrow(ME)
 end
end
end

function write_status(folder,status)
f=fopen(fullfile(folder,'pipeline_status.json'),'w');assert(f>=0);
c=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(status,PrettyPrint=true));
end
