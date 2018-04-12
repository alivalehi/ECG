function [ X ] = pi2twopi( X )
%PI2TWOPI Summary of this function goes here
%   Detailed explanation goes here
if(X<0)
    X = X+2*pi;
end
end

