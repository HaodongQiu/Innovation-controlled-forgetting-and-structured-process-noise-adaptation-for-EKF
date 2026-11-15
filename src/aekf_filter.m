function [x_est,s2nr,noise_var,history,metrics,diagnostic] = aekf_filter(M,K_undamaged,d,v,a,F,dt,noise_seed,noise_std,cfg,method,truth_k)
%AEKF_FILTER Joint response/stiffness EKF for the six manuscript configurations.
% The covariance adaptation, explicit Newmark transition and numerical
% safeguards are retained from the implementation used in the paper.
% Histories/metrics are passive; plotting is implemented separately.
% cfg.response_alpha selects 0.80 for fixed configurations, 0.95 otherwise.
dof=size(M,1); Nt=dof; Nx=4*dof; N=size(a,2); N0=dof;
beta=0; gamma=0.5; I_dof=eye(dof);
qa_floor_factor=0.01;
measured_acceleration_dofs=1:dof;
story_K_undamaged=extract_story_stiffness(K_undamaged);
noise_stream=RandStream('mt19937ar','Seed',noise_seed);
z_obs=zeros(dof,N); noise_var=zeros(dof,1); s2nr=zeros(dof,1);
for j=1:dof
    noise=noise_std*randn(noise_stream,1,N);
    z_obs(j,:)=a(j,:)+noise;
    noise_var(j)=var(noise,0,2);
    s2nr(j)=20*log10(norm(a(j,:))/norm(noise));
end
% Continuous self-tuning architecture:
% 1) the EKF innovation covariance SK is used only for the Kalman update;
% 2) an empirical innovation covariance supplies a meaningful NIS detector;
% 3) global change severity schedules stiffness covariance memory.
% No direct variance feedback, covariance jump, or localization is used.
alpha_Q_state       = cfg.response_alpha; % Fixed factor for displacement/velocity/acceleration
alpha_Q_initial     = 0.8; % Initial factor for every stiffness
alpha_Q_nominal     = 0.95; % Healthy-operation stiffness factor
alpha_Q_fast        = 0.8; % Smallest stiffness factor during a structural change
innovation_fast_rho = 0.98;    % Approximately 50-sample covariance memory
innovation_baseline_rho = 0.9995; % Slow healthy-baseline covariance memory
innovation_warmup_rho = 0.99;  % Faster baseline learning during initialization
innovation_warmup_time = 5.0;  % Seconds before damage adaptation is enabled
covariance_change_threshold = 1.15; % Fast/baseline covariance ratio threshold
covariance_change_scale = 1.50; % Additional ratio giving full severity
alpha_attack_rate   = 0.20;    % Rate of movement toward fast adaptation
alpha_recovery_rate = 0.002;   % Slow return toward the nominal factor

if ~(0 < alpha_Q_fast && alpha_Q_fast <= alpha_Q_initial && ...
        alpha_Q_initial < 1 && 0 < alpha_Q_state && alpha_Q_state < 1)
    error('Require 0 < alpha_Q_fast <= alpha_Q_initial < 1.');
end
if innovation_fast_rho <= 0 || innovation_fast_rho >= 1 || ...
        innovation_baseline_rho <= innovation_fast_rho || ...
        innovation_baseline_rho >= 1 || ...
        innovation_warmup_rho <= 0 || innovation_warmup_rho >= 1 || ...
        innovation_warmup_time <= 0 || ...
        covariance_change_threshold <= 1 || covariance_change_scale <= 0 || ...
        alpha_attack_rate <= 0 || alpha_attack_rate > 1 || ...
        alpha_recovery_rate <= 0 || alpha_recovery_rate > 1
    error('Invalid continuous self-tuning parameter.');
end

assert(ismember(method.controller,{'fixed','global'}));
if strcmp(method.controller,'fixed'),alpha_Q_initial=method.fixed_alpha;end

% Prescribed initialization. Response scales are fixed constants,
% independent of the current case response and observation layout.
noise_var_all_dofs = zeros(dof,1);
noise_var_all_dofs(measured_acceleration_dofs) = noise_var;
% Fixed prescribed response scales; reused unchanged for every case.
% Fixed prescribed values are supplied in data/models; see docs/SETTINGS.md.
response_acc_variance=cfg.response_acc_variance(:);
assert(numel(response_acc_variance)==dof && all(response_acc_variance>0),'Invalid prescribed prior');

q_acceleration=10*response_acc_variance;
Ba=[beta*dt^2*I_dof;gamma*dt*I_dof;I_dof];

% A 1% per sqrt(second) stiffness diffusion scale is used initially.
% below changes it multiplicatively; the limits are safety bounds, not
% additive covariance injections.
q_theta_nominal = (0.01*story_K_undamaged).^2 * dt;
q_theta_min     = (0.001*story_K_undamaged).^2 * dt;
q_theta_max     = (0.05*story_K_undamaged).^2 * dt;

% The initial stiffness estimate is deliberately 50% high, so its prior
% standard deviation is also set to 50% of nominal stiffness.
p_acceleration = response_acc_variance;
p_velocity = max((10*dt)^2*p_acceleration, 1e-6*ones(dof,1));
p_displacement = max((0.5*(10*dt)^2)^2*p_acceleration, ...
    1e-10*ones(dof,1));
p_theta = (0.50*story_K_undamaged).^2;

p = [
    p_displacement', ...
    p_velocity', ...
    p_acceleration', ...
    p_theta'
];

% Measurement covariance for measured channels only
R_base = diag(noise_var);

%% ========================================================================
%  10. Mass matrix operation
%  ========================================================================

% More stable than inv(M)
Minv = M \ eye(dof);

%% ========================================================================
%  11. Constant Rayleigh damping matrix
%  ========================================================================

K0 = K_undamaged;

omega0 = compute_undamped_eig(M, K0);

alpha_R0 = omega0(1)/omega0(2);

beta_R0 = ...
    (1-alpha_R0)*0.01 / ...
    (omega0(2)-alpha_R0*omega0(1));

a_Rayleigh0 = 2*omega0(1)*omega0(2)*beta_R0;
b_Rayleigh0 = 2*beta_R0;

D_const = a_Rayleigh0*M + b_Rayleigh0*K0;



Q=blkdiag(Ba*diag(q_acceleration)*Ba',diag(q_theta_nominal));
R=R_base; P=diag(p);
x_est=[zeros(3*dof,1);1.5*story_K_undamaged];
alpha_Q_story_current=alpha_Q_initial*ones(Nt,1);
q_theta_current=q_theta_nominal;
innovation_covariance_fast=R_base; innovation_covariance_baseline=R_base;
innovation_covariance_initialized=false;
innovation_warmup_samples=min(N,max(N0+1,round(innovation_warmup_time/dt)));
numerical_limit=1e100;
story_difference=sparse(eye(dof)-diag(ones(dof-1,1),-1));
observed_indices = 2*dof + measured_acceleration_dofs;
theta_indices = 3*dof + (1:Nt);
% These operators are constant throughout filtering.
Minv_story_transpose = sparse(Minv*story_difference');
Minv_D = sparse(Minv*D_const);
Minv_F = Minv*F;
I_sparse = speye(dof);
stiffness_rows=[(1:dof)';(1:dof-1)';(2:dof)'];
stiffness_cols=[(1:dof)';(2:dof)';(1:dof-1)'];
lower_options=struct('LT',true);
upper_options=struct('UT',true);
newmark_noise_weights = [beta*dt^2; gamma*dt; 1];
response_q_indices = zeros(dof,9);
response_q_coefficients = zeros(1,9);
for row_block=1:3
    for col_block=1:3
        block_id=(row_block-1)*3+col_block;
        response_q_indices(:,block_id)=sub2ind([Nx,Nx], ...
            (row_block-1)*dof+(1:dof)',(col_block-1)*dof+(1:dof)');
        response_q_coefficients(block_id)=newmark_noise_weights(row_block)*newmark_noise_weights(col_block);
    end
end
theta_q_indices=sub2ind([Nx,Nx],theta_indices,theta_indices);
qa_current=q_acceleration;

% Small matrices run faster without repeated BLAS thread scheduling.
% Restore the caller's thread setting immediately after the timed loop.
previous_thread_count = maxNumCompThreads(1);
thread_cleanup = onCleanup(@() maxNumCompThreads(previous_thread_count));

% BEGIN OBSERVATION-ONLY INITIALIZATION
sample_indices=double(unique([1:cfg.history_stride:N,N,cfg.half_N,cfg.half_N+1]));
sample_indices=sample_indices(sample_indices>=1 & sample_indices<=N);
sample_slot=zeros(1,N);sample_slot(sample_indices)=1:numel(sample_indices);
history.time=(sample_indices-1)*dt;
history.theta=nan(Nt,numel(sample_indices));history.alpha=history.theta;history.q=history.theta;
history.theta_sd=nan(Nt,numel(sample_indices));
history.response_estimate=nan(3*dof,numel(sample_indices));
history.response_true=history.response_estimate;
history.ratio=nan(1,numel(sample_indices));history.severity=history.ratio;history.nis=history.ratio;history.empirical_normalized_nis=history.ratio;
metrics.error_energy=zeros(3*dof,size(cfg.windows,1));metrics.reference_energy=metrics.error_energy;
metrics.count=zeros(1,size(cfg.windows,1));
metrics.acceleration_noise_energy=zeros(dof,size(cfg.windows,1));
for wi=1:size(cfg.windows,1)
    take=cfg.windows(wi,1):cfg.windows(wi,2);
    metrics.acceleration_noise_energy(:,wi)=sum((z_obs(:,take)-a(:,take)).^2,2);
end
metrics.stiffness_error_energy=zeros(dof,size(cfg.windows,1));
metrics.damage_error_energy=metrics.stiffness_error_energy;
diagnostic=struct('complete',false,'failure_message','','last_sample',0,'min_stiffness',Inf, ...
    'minimum_scaled_Q_eigenvalue',Inf,'initial_alpha',alpha_Q_story_current, ...
    'initial_P_diagonal',diag(P),'initial_Q_diagonal',diag(Q), ...
    'initial_state',x_est,'R_diagonal',diag(R),'positive_definiteness_checks',0);
Qxx_current=Q(1:3*dof,1:3*dof);
last_violation=cfg.half_N;last_forced_violation=cfg.half_N;
last_damage_band_violation=cfg.half_N;last_forced_damage_band_violation=cfg.half_N;
diagnostic.initialization_uses_future_response=false;
diagnostic.controller=method.controller;
diagnostic.response_forgetting_factor=alpha_Q_state;
diagnostic.fixed_stiffness_forgetting_factor=method.fixed_alpha;
diagnostic.nis_controls_forgetting=~strcmp(method.controller,'fixed');
diagnostic.response_mapping=method.response_mapping;
diagnostic.parameter_independence=method.parameter_independence;
history.core_step_seconds=nan(1,N);
history.full_empirical_normalized_nis=nan(1,N);
history.full_applied_alpha=nan(1,N);
history.full_severity=nan(1,N);
benchmark_timer=tic;
try
% END OBSERVATION-ONLY INITIALIZATION
for data_index=1:N
step_timer=tic;
Q_previous=Q; empirical_nis=NaN;
story_k = x_est(theta_indices);
K = sparse(stiffness_rows,stiffness_cols, ...
    [story_k+[story_k(2:end);0]; -story_k(2:end); -story_k(2:end)],dof,dof);

% Constant damping matrix
D = D_const;

%% ---------------------------------------------------------------
%  Nonlinear state-transition function
%  ---------------------------------------------------------------

Minv_K = sparse(Minv*K);
acceleration_d = -Minv_K;
acceleration_v = -(Minv_D + dt*Minv_K);
acceleration_a = -(dt^2*(0.5-beta)*Minv_K + dt*(1-gamma)*Minv_D);
A = [I_sparse + beta*dt^2*acceleration_d, ...
     dt*I_sparse + beta*dt^2*acceleration_v, ...
     dt^2*(0.5-beta)*I_sparse + beta*dt^2*acceleration_a;
     gamma*dt*acceleration_d, ...
     I_sparse + gamma*dt*acceleration_v, ...
     dt*(1-gamma)*I_sparse + gamma*dt*acceleration_a;
     acceleration_d, acceleration_v, acceleration_a];

% The current estimate is for data sample data_index, so this
% uses the Newmark new-time load F_{n+1}, not the previous load F_n.
B = [
    beta*dt^2*Minv_F(:, data_index);
    gamma*dt*Minv_F(:, data_index);
    Minv_F(:, data_index)
];



%% ---------------------------------------------------------------
%  Jacobian matrix
%  ---------------------------------------------------------------

J = [
    A, sparse(3*dof, Nt);
    sparse(Nt, 3*dof), speye(Nt)
];

x_previous=x_est(1:3*dof);
dp_previous=x_previous(1:dof)+dt*x_previous(dof+1:2*dof)+dt^2*(0.5-beta)*x_previous(2*dof+1:3*dof);
Gtheta=-Minv_story_transpose*sparse(1:dof,1:dof,story_difference*dp_previous,dof,dof);
J(1:3*dof,3*dof+(1:Nt))=[beta*dt^2*Gtheta;gamma*dt*Gtheta;Gtheta];

%% ---------------------------------------------------------------
%  Prediction
%  ---------------------------------------------------------------
f = [
    A*x_est(1:3*dof) + B;
    x_est(3*dof+1:Nx)
];
u_pre = f;

Pp = J*P*J' + Q;

% Remove small numerical asymmetry
Pp = 0.5*(Pp + Pp');

% Stop this model before non-finite Pp contaminates SK.
matrix_peak = max(abs(Pp),[],'all','includenan');
if ~isfinite(matrix_peak) || matrix_peak > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Predicted covariance Pp contains NaN/Inf or overflow.');
end

%% ---------------------------------------------------------------
%  Measurement update
%  ---------------------------------------------------------------

Innovation = z_obs(:, data_index) - u_pre(observed_indices);

if any(~isfinite(Innovation)) || max(abs(Innovation)) > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Innovation contains NaN/Inf or overflow.');
end

SK = Pp(observed_indices,observed_indices) + R;

% Remove small numerical asymmetry
SK = 0.5*(SK + SK');

matrix_peak = max(abs(SK),[],'all','includenan');
if ~isfinite(matrix_peak) || matrix_peak > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Innovation covariance SK contains NaN/Inf or overflow.');
end

% A non-positive-definite SK indicates that this model's covariance
% has become numerically invalid. Stop this model only.
[L_SK, chol_flag] = chol(SK, 'lower');

if chol_flag ~= 0
    error('AEKF:ModelDivergence', ...
        'Innovation covariance SK is not positive definite.');
end

% Kalman gain
whitened_cross = linsolve(L_SK,Pp(:,observed_indices)',lower_options)';
Kg = linsolve(L_SK',whitened_cross',upper_options)';

matrix_peak = max(abs(Kg),[],'all','includenan');
if ~isfinite(matrix_peak) || matrix_peak > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Kalman gain contains NaN/Inf or overflow.');
end


%% Global or fixed forgetting; no localization.
if strcmp(method.controller,'fixed')
    normalized_innovation=linsolve(L_SK,Innovation,lower_options);
    nis_current=normalized_innovation'*normalized_innovation;
    alpha_Q_story_current(:)=method.fixed_alpha;covariance_change_ratio=NaN;change_severity=0;
else
%  Empirical-NIS detector for global dynamic forgetting
%  ---------------------------------------------------------------
% SK contains deliberately conservative state/process uncertainty,
% so its NIS is retained only as a filter-consistency diagnostic.
normalized_innovation = linsolve(L_SK,Innovation,lower_options);
nis_current = normalized_innovation' * normalized_innovation;

innovation_outer = Innovation*Innovation';
if ~innovation_covariance_initialized
    % Online bootstrap: diagonal initialization is positive definite
    % and uses no future samples.
    bootstrap_variance = max(Innovation.^2, diag(R_base));
    innovation_covariance_baseline = ...
        diag(bootstrap_variance) + R_base;
    innovation_covariance_fast = innovation_covariance_baseline;
    innovation_covariance_initialized = true;
end

innovation_covariance_baseline = 0.5 * ...
    (innovation_covariance_baseline + ...
     innovation_covariance_baseline');
L_empirical = chol(innovation_covariance_baseline, 'lower');
empirical_normalized_innovation = linsolve(L_empirical,Innovation,lower_options);
empirical_nis = empirical_normalized_innovation' * ...
    empirical_normalized_innovation;

innovation_covariance_fast = ...
    innovation_fast_rho*innovation_covariance_fast + ...
    (1-innovation_fast_rho)*innovation_outer;
innovation_covariance_fast = 0.5 * ...
    (innovation_covariance_fast + innovation_covariance_fast');

covariance_change_ratio = trace( ...
    innovation_covariance_baseline \ ...
    innovation_covariance_fast)/N0;
covariance_change_ratio = max(real(covariance_change_ratio), 0);

if data_index <= innovation_warmup_samples
    change_severity = 0;
    baseline_rho_current = innovation_warmup_rho;
else
    change_excess = max(covariance_change_ratio - ...
        covariance_change_threshold, 0);
    change_severity = min(change_excess / ...
        covariance_change_scale, 1);

    % The healthy baseline continues to learn, but its update slows
    % continuously as evidence of a structural change increases.
    baseline_learning_fraction = ...
        (1-innovation_baseline_rho)*(1-change_severity);
    baseline_rho_current = 1-baseline_learning_fraction;
end

innovation_covariance_baseline = ...
    baseline_rho_current*innovation_covariance_baseline + ...
    (1-baseline_rho_current)*innovation_outer;
innovation_covariance_baseline = 0.5 * ...
    (innovation_covariance_baseline + ...
     innovation_covariance_baseline');

% Global dynamic forgetting: every stiffness receives the SAME target.
% No parameter-sensitivity scores, spatial weights or localization are used.
alpha_Q_story_target = alpha_Q_nominal-change_severity*(alpha_Q_nominal-alpha_Q_fast)*ones(Nt,1);

for parameter_id = 1:Nt
    if alpha_Q_story_target(parameter_id) < ...
            alpha_Q_story_current(parameter_id)
        alpha_rate = alpha_attack_rate;
    else
        alpha_rate = alpha_recovery_rate;
    end
    alpha_Q_story_current(parameter_id) = ...
        alpha_Q_story_current(parameter_id) + alpha_rate * ...
        (alpha_Q_story_target(parameter_id) - ...
         alpha_Q_story_current(parameter_id));
end
alpha_Q_story_current = min(max(alpha_Q_story_current, ...
    alpha_Q_fast), alpha_Q_nominal);

end

%% ---------------------------------------------------------------
%  Continuous multiplicative process-noise update
%  ---------------------------------------------------------------
% Only the diagonal of this rank-one matrix is used below.
state_correction = Kg*Innovation;
qdata = state_correction.^2;

% PSD-preserving covariance matching for the complete augmented Q.
alpha_Q_vector = [
    alpha_Q_state*ones(3*dof,1);
    alpha_Q_story_current
];
% Independent acceleration disturbances propagate through Newmark:
% [delta d; delta v; delta a] = Ba * delta a. Stiffness diffusion
% is a separate driving process. Retain all posterior P correlations.
qa_current=alpha_Q_state*qa_current+(1-alpha_Q_state)*...
    max(qdata(2*dof+(1:dof)),qa_floor_factor*noise_var_all_dofs);
% Geometric covariance matching with one common stiffness forgetting factor.
% Final variance bounds are applied below.
q_theta_data = max(qdata(theta_indices), q_theta_min);
log_q_blended = alpha_Q_story_current.*log(q_theta_current) + ...
    (1-alpha_Q_story_current).*log(q_theta_data);
q_theta_current = exp(log_q_blended);
q_theta_current = min(max(q_theta_current, q_theta_min), ...
    q_theta_max);

% Same final block-diagonal Q as the original congruence scaling.
Q(response_q_indices) = response_q_coefficients .* qa_current;
Q(theta_q_indices) = q_theta_current;
if method.parameter_independence
    % M2: parameter driving disturbances independent; Q_x_theta=0.
    if strcmp(method.response_mapping,'unmapped')
        cx=state_correction(1:3*dof);
        Qxx_current=alpha_Q_state*Qxx_current+(1-alpha_Q_state)*(cx*cx');
        ai=2*dof+(1:dof);
        missing_variance=max(qa_current-diag(Qxx_current(ai,ai)),0);
        Qxx_current(ai,ai)=Qxx_current(ai,ai)+diag(missing_variance);
        Q(1:3*dof,1:3*dof)=Qxx_current;
    end
else

% Without M2, retain correlated parameter driving disturbances.
% All diagonal parameter variances use the same bounded geometric law.
memory=sqrt(alpha_Q_vector);fresh=sqrt(1-alpha_Q_vector).*state_correction;
raw_Q=(memory*memory').*Q_previous+fresh*fresh';
ai=2*dof+(1:dof);ti=theta_indices;
extra=max(qa_current-diag(raw_Q(ai,ai)),0);
raw_Q(ai,ai)=raw_Q(ai,ai)+diag(extra);
theta_scale=sqrt(q_theta_current./max(diag(raw_Q(ti,ti)),realmin));
scale_vector=[ones(3*dof,1);theta_scale];
raw_Q=(scale_vector*scale_vector').*raw_Q;
if strcmp(method.response_mapping,'mapped')
    % The state mapping and cross-covariance must be consistent for PSD.
    % This is the symmetric-whitening construction in the supplied 20-DOF
    % response-only source; parameter/parameter covariance is retained.
    Qa=0.5*(raw_Q(ai,ai)+raw_Q(ai,ai)');qa=diag(Qa);s=max(qa);
    [V,Lambda]=eig(Qa/s);lambda=real(diag(Lambda));
    assert(all(lambda>0),'Acceleration principal covariance is not positive definite');
    W=diag(sqrt(qa/s))*V*diag(1./sqrt(lambda))*V';
    qat=W*raw_Q(ai,ti);
    Q=zeros(Nx);Q(1:3*dof,1:3*dof)=Ba*diag(qa)*Ba';
    Q(1:3*dof,ti)=Ba*qat;Q(ti,1:3*dof)=Q(1:3*dof,ti)';Q(ti,ti)=raw_Q(ti,ti);
else
    Q=raw_Q;
end
Q=0.5*(Q+Q');

end


matrix_peak = max(abs(Q),[],'all','includenan');
if ~isfinite(matrix_peak) || matrix_peak > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Adaptive process covariance Q contains NaN/Inf or overflow.');
end

%% ---------------------------------------------------------------
%  State update
%  ---------------------------------------------------------------

x_est = ...
    u_pre + state_correction;

if any(~isfinite(x_est)) || ...
        max(abs(x_est)) > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Updated state contains NaN/Inf or overflow.');
end

%% ---------------------------------------------------------------
%  Joseph-form covariance update
%  ---------------------------------------------------------------

% Complete-square form of the SAME Joseph covariance update:
% C=Pp*H', W=C/L', SK=L*L', E=Kg*L-W.
% (I-Kg*H)*Pp*(I-Kg*H)' + Kg*R*Kg' = Pp-W*W'+E*E'.
% E retains the gain-solve roundoff correction rather than assuming zero.
gain_solve_residual = Kg*L_SK-whitened_cross;
P = Pp - whitened_cross*whitened_cross' + ...
    gain_solve_residual*gain_solve_residual';

% Remove numerical asymmetry
P = ...
    0.5*(P + P');

matrix_peak = max(abs(P),[],'all','includenan');
if ~isfinite(matrix_peak) || matrix_peak > numerical_limit
    error('AEKF:ModelDivergence', ...
        'Posterior covariance P contains NaN/Inf or overflow.');
end


history.core_step_seconds(data_index)=toc(step_timer);
history.full_empirical_normalized_nis(data_index)=empirical_nis/N0;
history.full_applied_alpha(data_index)=alpha_Q_story_current(1);
history.full_severity(data_index)=change_severity;
% BEGIN OBSERVATION-ONLY RECORDING
diagnostic.last_sample=data_index;
diagnostic.min_stiffness=min(diagnostic.min_stiffness,min(x_est(theta_indices)));
active=find(data_index>=cfg.windows(:,1) & data_index<=cfg.windows(:,2));
reference_sample=[d(:,data_index);v(:,data_index);a(:,data_index)];
err=x_est(1:3*dof)-reference_sample;
metrics.error_energy(:,active)=bsxfun(@plus,metrics.error_energy(:,active),err.^2);
metrics.reference_energy(:,active)=bsxfun(@plus,metrics.reference_energy(:,active),reference_sample.^2);
metrics.count(active)=metrics.count(active)+1;
if data_index>cfg.half_N,kt=truth_k;else,kt=story_K_undamaged;end
ke=100*(x_est(theta_indices)./kt-1);
de=100*(x_est(theta_indices)-kt)./story_K_undamaged;
metrics.stiffness_error_energy(:,active)=bsxfun(@plus,metrics.stiffness_error_energy(:,active),ke.^2);
metrics.damage_error_energy(:,active)=bsxfun(@plus,metrics.damage_error_energy(:,active),de.^2);
if data_index>cfg.half_N && max(abs(de))>cfg.damage_band_pp
    last_damage_band_violation=data_index;
    if data_index<=cfg.forced_samples,last_forced_damage_band_violation=data_index;end
end
if data_index==cfg.forced_samples
    metrics.forced_end_theta=x_est(theta_indices);
end
if data_index>cfg.half_N
    relative_error=max(abs(x_est(theta_indices)./truth_k-1))*100;
    if relative_error>5
        last_violation=data_index;
        if data_index<=cfg.forced_samples,last_forced_violation=data_index;end
    end
end
slot=sample_slot(data_index);
if slot>0
    history.theta(:,slot)=x_est(theta_indices);history.alpha(:,slot)=alpha_Q_story_current;
    history.q(:,slot)=q_theta_current;
    history.theta_sd(:,slot)=sqrt(max(diag(P(theta_indices,theta_indices)),0));
    history.response_estimate(:,slot)=x_est(1:3*dof);
    history.response_true(:,slot)=reference_sample;
    history.ratio(slot)=covariance_change_ratio;history.severity(slot)=change_severity;history.nis(slot)=nis_current;history.empirical_normalized_nis(slot)=empirical_nis/N0;
end
if mod(data_index,2000)==0 || data_index==N
    qdiag=diag(Q);keep=qdiag>0;
    scale=sqrt(qdiag(keep));normalized_Q=Q(keep,keep)./(scale*scale');
    smallest=min(real(eig(0.5*(normalized_Q+normalized_Q'))));
    diagnostic.minimum_scaled_Q_eigenvalue=min(diagnostic.minimum_scaled_Q_eigenvalue,smallest);
    diagnostic.positive_definiteness_checks=diagnostic.positive_definiteness_checks+1;
    assert(smallest>-1e-7,'Comparison covariance lost positive semidefiniteness');
    assert(max(abs(diag(Q(3*dof+1:4*dof,3*dof+1:4*dof))./q_theta_current-1))<1e-10,'Parameter marginal variance law changed');
end
% END OBSERVATION-ONLY RECORDING
end
diagnostic.complete=true;
catch ME
    diagnostic.failure_message=getReport(ME,'extended','hyperlinks','off');
end
diagnostic.seconds=toc(benchmark_timer);
diagnostic.last_violation_index=last_violation;
diagnostic.last_forced_violation_index=last_forced_violation;
diagnostic.final_Q=Q;
diagnostic.last_damage_band_violation=last_damage_band_violation;
diagnostic.last_forced_damage_band_violation=last_forced_damage_band_violation;
metrics.stiffness_window_rmse_pct=bsxfun(@rdivide,sqrt(metrics.stiffness_error_energy),sqrt(metrics.count));
metrics.damage_window_rmse_pp=bsxfun(@rdivide,sqrt(metrics.damage_error_energy),sqrt(metrics.count));
metrics.rmse=bsxfun(@rdivide,sqrt(metrics.error_energy),sqrt(metrics.count));
metrics.reference_rms=bsxfun(@rdivide,sqrt(metrics.reference_energy),sqrt(metrics.count));
metrics.normalized_rmse_pct=100*sqrt(metrics.error_energy./metrics.reference_energy);
metrics.normalized_rmse_pct(metrics.reference_energy==0)=NaN;
clear thread_cleanup;
end

function k_story = extract_story_stiffness(K)

    dof = size(K,1);
    k_story = zeros(dof,1);
    k_story(2:dof) = -diag(K,1);
    k_story(1) = K(1,1)-k_story(2);

end

function omega = compute_undamped_eig(M, K)

    eigval = real(eig(K, M));

    % Retain the positive structural vibration modes.
    eig_scale = max(abs(eigval));
    eig_tol = max(1, eig_scale)*1e-10;
    eigval = sort(eigval(eigval > eig_tol));

    if numel(eigval) < 2
        error('At least two positive vibration modes are required.');
    end

    % Return circular frequencies in rad/s for the Rayleigh formulas.
    omega = sqrt(eigval);

end
