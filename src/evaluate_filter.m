function result = evaluate_filter(model,spec,cfg,method,generator,d,v,a,F,truth_info,signature)
%EVALUATE_FILTER Execute one configuration and apply paper success criterion.
result=struct('protocol','MSSP_public_v1','signature',signature, ...
 'model_id',model.id,'case_id',spec.case_id,'method',method, ...
 'generator',generator,'cfg',cfg,'spec',spec,'nominal_k',model.k0, ...
 'true_k',spec.true_k,'completed',false,'success',false,'failure_message','');
timer=tic;
try
 [xf,result.snr_db,result.noise_variance,result.trace,result.metrics,result.diagnostic]= ...
  aekf_filter(model.M,model.K0,d,v,a,F,model.dt,spec.noise_seed,cfg.noise_std,cfg,method,spec.true_k);
 result.final_state=xf;result.final_true_state=[d(:,end);v(:,end);a(:,end)];
 result.completed=true; % A returned divergence is a finalized comparison.
 result.success=result.diagnostic.complete && all(isfinite(xf)) && all(xf(3*model.n+1:end)>0);
 result.failure_message=result.diagnostic.failure_message;
 if result.diagnostic.complete && ~result.success
  result.failure_message='Nonpositive final stiffness: parameter divergence.';
 end
catch ME
 % Programming, input or filesystem errors must stop the campaign.
 % Numerical divergence handled inside aekf_filter returns normally above.
 rethrow(ME)
end
result.truth_info=truth_info;result.total_seconds=toc(timer);
end
