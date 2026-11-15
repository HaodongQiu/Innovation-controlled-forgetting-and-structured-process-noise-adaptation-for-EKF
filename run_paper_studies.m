function run_paper_studies(workers)
%RUN_PAPER_STUDIES All numerical campaigns used in the manuscript.
% Plotting is a separate step: run plot_paper_results afterwards.
if nargin<1,workers=4;end
run_monte_carlo(workers);
run_noise_sensitivity(workers);
run_response_generator_study(workers);
end
