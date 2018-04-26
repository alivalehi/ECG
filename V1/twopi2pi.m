function [ X ] = twopi2pi( X )
%PI2TWOPI Summary of this function goes here
%   Detailed explanation goes here
for i =1:length(X)
if(X>pi)
    X = (X-2*pi);
end
end

end

