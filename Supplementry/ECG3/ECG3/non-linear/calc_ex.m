function [snrdbc, Ex2] = calc_ex(mu, sig, pai)
%to present correct snr
[p, t] = size(mu);   Ex2 = 0;
for i = 1: t,     Ex2 = Ex2 + pai(i)*(mu(:,i).^2 + diag(sig(:,:,i))); end
snrc = sum(Ex2)/p; Ex2 = snrc; snrdbc = 10*log10(snrc);
return
