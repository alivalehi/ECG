function [ X ] = twopi2pi( X )
%PI2TWOPI Summary of this function goes here
%   Detailed explanation goes here
if(X>pi)
    X = (X-2*pi);
end
end

