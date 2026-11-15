function run_response_generator_study(workers)
%RUN_RESPONSE_GENERATOR_STUDY Proposed filter; Newmark, ode45 and exponential.
if nargin<1,workers=4;end
for generator={'newmark','ode45','expm_foh'}
 run_monte_carlo(workers,{'Pirelli39'},.01,generator{1},1:100,{'M12-dynamic'});
end
end
