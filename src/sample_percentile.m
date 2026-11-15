function q = sample_percentile(x,p)
%SAMPLE_PERCENTILE MATLAB-compatible midpoint sample percentile; no toolbox.
x=sort(x(isfinite(x)));x=x(:).';n=numel(x);
if n==0,q=nan(size(p));return;end
position=min(max(n*p/100+.5,1),n);lo=floor(position);hi=ceil(position);
q=x(lo)+(position-lo).*(x(hi)-x(lo));
q=reshape(q,size(p));
end
