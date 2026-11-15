function [d,v,a,F,info]=continuous_response(model,damaged_k,generator,cfg)
% Continuous linear structural dynamics, with identical FOH input in both solvers.
% At the switch, d/v are continuous. At the exact boundary the saved a uses
% the healthy (left) limit, preserving the archived sample convention.
n=model.n;N=model.N;h=model.dt;split=model.half_N;
E=eye(n)-diag(ones(n-1,1),-1);Kd=E'*diag(damaged_k)*E;
A0=[zeros(n) eye(n);-(model.M\model.K0) -(model.M\model.D)];
A1=[zeros(n) eye(n);-(model.M\Kd) -(model.M\model.D)];
B=[zeros(n,1);-ones(n,1)];u=model.ag_full;t=(0:N-1)*h;
old=maxNumCompThreads(1);cleanup=onCleanup(@() maxNumCompThreads(old));timer=tic;
switch generator
    case 'expm_foh'
        [P0,G0,H0]=foh_operators(A0,B,h);[P1,G1,H1]=foh_operators(A1,B,h);
        z=zeros(2*n,N);
        for j=1:N-1
            if j<split,P=P0;G=G0;H=H0;else,P=P1;G=G1;H=H1;end
            z(:,j+1)=P*z(:,j)+G*u(j)+H*((u(j+1)-u(j))/h);
        end
    case 'ode45'
        opts=odeset('RelTol',cfg.ode_rel_tol, ...
            'AbsTol',[cfg.ode_abs_d*ones(n,1);cfg.ode_abs_v*ones(n,1)], ...
            'MaxStep',cfg.ode_max_step);
        [~,zh]=ode45(@(tt,zz) rhs(tt,zz,A0,B,u,h,N),t(1:split),zeros(2*n,1),opts);
        [~,zd]=ode45(@(tt,zz) rhs(tt,zz,A1,B,u,h,N),t(split:N),zh(end,:).',opts);
        z=[zh;zd(2:end,:)].';
    otherwise
        error('Unknown truth generator: %s',generator);
end
d=z(1:n,:);v=z(n+1:2*n,:);F=-(model.M*ones(n,1))*u;
a=zeros(n,N);
a(:,1:split)=model.M\(F(:,1:split)-model.D*v(:,1:split)-model.K0*d(:,1:split));
a(:,split+1:N)=model.M\(F(:,split+1:N)-model.D*v(:,split+1:N)-Kd*d(:,split+1:N));
assert(all(isfinite(z),'all')&&all(isfinite(a),'all'),'Truth generator produced nonfinite values');
rd=diff(d,1,2)-h*v(:,1:end-1)-.5*h^2*a(:,1:end-1);
rv=diff(v,1,2)-.5*h*(a(:,1:end-1)+a(:,2:end));
% Exclude the switch interval for a smooth-segment kinematic diagnostic.
keep=true(1,N-1);keep(split)=false;
info=struct('generator',generator,'seconds',toc(timer),'input_interpolation','piecewise_linear_FOH', ...
    'damage_switch_time',model.damage_switch_time,'first_damaged_sample_time',model.first_damaged_sample_time, ...
    'newmark_displacement_defect_rms_m',sqrt(mean(rd(:,keep).^2,2)), ...
    'newmark_velocity_defect_rms_mps',sqrt(mean(rv(:,keep).^2,2)), ...
    'defect_excludes_switch_interval',true,'final_true_state',[d(:,end);v(:,end);a(:,end)]);
info.relative_equilibrium_residual=norm(model.M*a+model.D*v+ ...
    [model.K0*d(:,1:split),Kd*d(:,split+1:end)]-F,'fro')/max(norm(F,'fro'),realmin);
end

function [P,G,H]=foh_operators(A,B,h)
% Augmented state [z; input; input slope]; exact for linear input on each step.
p=size(A,1);L=zeros(p+2);L(1:p,1:p)=A;L(1:p,p+1)=B;L(p+1,p+2)=1;
X=expm(h*L);P=X(1:p,1:p);G=X(1:p,p+1);H=X(1:p,p+2);
end

function dz=rhs(t,z,A,B,u,h,N)
q=min(max(t/h,0),N-1);j=min(floor(q)+1,N-1);f=q-(j-1);
ug=(1-f)*u(j)+f*u(j+1);dz=A*z+B*ug;
end
