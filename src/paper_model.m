function model = paper_model(model_id)
%PAPER_MODEL Load SI model matrices and the supplied earthquake acceleration.
arguments
 model_id (1,:) char {mustBeMember(model_id,{'Pirelli39','ASCE4'})}
end
root=fileparts(fileparts(mfilename('fullpath')));
model=load(fullfile(root,'data','models',[model_id '_model.mat']));
model.id=model_id;
for name={'n','N','half_N','forced_samples','dt','zeta'}
 model.(name{1})=double(model.(name{1}));
end
model.k0=model.k0(:);model.response_acc_variance=model.response_acc_variance(:);
input=load(fullfile(root,'data','ElCentro_time_history_dt0001.mat'));
model.ag=input.acc_new(:).';
model.ag_full=[model.ag zeros(1,model.N-numel(model.ag))];
omega=sort(sqrt(real(eig(model.K0,model.M))));
b=2*model.zeta/(omega(1)+omega(2));
model.D=b*omega(1)*omega(2)*model.M+b*model.K0;
model.frequencies_Hz=omega/(2*pi);
model.damage_switch_time=(model.half_N-1)*model.dt;
model.first_damaged_sample_time=model.half_N*model.dt;
E=eye(model.n)-diag(ones(model.n-1,1),-1);
assert(norm(E'*diag(model.k0)*E-model.K0,'fro')<1e-10*norm(model.K0,'fro'));
assert(model.N==numel(model.ag)+50000 && model.half_N==floor(numel(model.ag)/2));
assert(max(abs(diff(input.time_new(:))-model.dt))<1e-10);
end
