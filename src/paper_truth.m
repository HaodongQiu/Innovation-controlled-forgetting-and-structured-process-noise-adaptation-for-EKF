function [d,v,a,F,info] = paper_truth(model,spec,generator,cfg)
%PAPER_TRUTH Select truth generation; the filter always uses Newmark.
if strcmp(generator,'newmark')
 timer=tic;[d,v,a,F]=newmark_response(model,spec.true_k);
 info=struct('generator','newmark','seconds',toc(timer));
else
 [d,v,a,F,info]=continuous_response(model,spec.true_k,generator,cfg);
end
end
