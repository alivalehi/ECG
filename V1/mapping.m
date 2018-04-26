function [ result ] = mapping(X,fard,closed,far,close,i,twopi_range)
try
if nargin < 7
    twopi_range = 0;
end
if(twopi_range)
    %Some code for adjusting range of azimuth which is -pi to pi
    X(end) = pi2twopi(X(end));
    fard = pi2twopi(fard);
    far = pi2twopi(far);
    closed = pi2twopi(closed);
    close = pi2twopi(close);
end
alpha1 = 100;%abs(far-close)/10 ;
if(fard>closed)
    result = logit([0,close,far,2*pi],[0,closed,fard,2*pi],X(i),alpha1);
else
     result = logit_reverse([0,close,far,2*pi],[2*pi,closed,fard,0],X(i),alpha1);
end    
    if(twopi_range)  
    result = twopi2pi(result);    
end
catch
   error('error in mapping function'); 
end
end

