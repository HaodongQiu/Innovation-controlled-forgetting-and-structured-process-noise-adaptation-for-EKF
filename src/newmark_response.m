function [d,v,a,F] = newmark_response(model,damaged_k)
%NEWMARK_RESPONSE Supplied explicit Newmark recurrence, beta=0, gamma=1/2.
% Relative structural responses and base-excitation load are in SI units.
n=model.n;N=model.N;dt=model.dt;
E=eye(n)-diag(ones(n-1,1),-1);Kd=E'*diag(damaged_k)*E;
F=zeros(n,N);F(:,1:numel(model.ag))=-(model.M*ones(n,1))*model.ag;
MinvF=model.M\F;MinvD=model.M\model.D;
MinvK0=model.M\model.K0;MinvKd=model.M\Kd;
d=zeros(n,N);v=d;a=d;a(:,1)=MinvF(:,1);
for j=1:N-1
 if j<model.half_N,MinvK=MinvK0;else,MinvK=MinvKd;end
 dp=d(:,j)+dt*v(:,j)+.5*dt^2*a(:,j);
 vp=v(:,j)+.5*dt*a(:,j);
 a(:,j+1)=MinvF(:,j+1)-MinvK*dp-MinvD*vp;
 d(:,j+1)=dp;v(:,j+1)=vp+.5*dt*a(:,j+1);
end
assert(all(isfinite(a),'all') && all(isfinite(v),'all') && all(isfinite(d),'all'), ...
 'Truth response diverged.');
end
