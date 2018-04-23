function [ X ] = pi2twopi( X )
%PI2TWOPI Summary of this function goes here
%   Detailed explanation goes here
for i =1:length(X)
if(X(i)<0)
    X(i) = X(i)+2*pi;
end
end
end

