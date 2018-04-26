function [ Y ] = predict( X,Craw,Ccooked,CloseIndex,FarIndex )

%PREDICT Summary of this function goes here
%   Input: X input cartesian
%          C contains matrix of centroids in spherical
%          CloseIndex  Close index
%          FarIndex
%  Output: Y transformed data
sph_all_2_ref  = X - repmat(Craw(1,:), size(X,1), 1);
[~,dimension] = size(Craw);
twopi_range=0;
for i = 2:dimension
    if(dimension==i)
        twopi_range =1;
    end
    X(:) = script_roation (X,mapping(X,Ccooked(FarIndex(i-1),i),Ccooked(CloseIndex(i-1),i),...
        Craw(FarIndex(i-1),i),Craw(CloseIndex(i-1),i),i-1,twopi_range),i);
end
Y = nsphere2cart(X);

