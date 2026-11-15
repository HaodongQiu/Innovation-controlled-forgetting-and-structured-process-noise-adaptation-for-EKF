function run_noise_sensitivity(workers)
%RUN_NOISE_SENSITIVITY Proposed Pirelli method, 100 matched cases x four noises.
if nargin<1,workers=4;end
for sigma=[.005 .01 .03 .05]
 run_monte_carlo(workers,{'Pirelli39'},sigma,'newmark',1:100,{'M12-dynamic'});
end
end
