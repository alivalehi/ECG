function [ T ] = Fisher( Y,delimiter)
%FISHER Summary of this function goes here
%[X11 X21 X31
% X12 X22 X32
% X13 X23 X33]
Sksigma = @(X) sum((X-repmat(mean(X,2),1,length(X)))*(X-repmat(mean(X,2),1,length(X)))',2);%the X which is passed has to be data correspond o one cluster
Nclusters = length(delimiter)-1;%0,500,1000,1500 => 3 clusters
for j=1:Nclusters
    input = Y(:,delimiter(j):delimiter(j+1));
    Sk(:,j) = Sksigma(input);
end
Nk = @(X) length(X);
mk = @(X) mean(X,2);
Sw = sum(Sk);
for i=1:Nclusters
    input = Y(:,delimiter(i):delimiter(i+1));
    Sb(:,i) = sum(Nk(input) * (mk(input)-mean(Y,2)) * (mk(input)-mean(Y,2))',2);
end
Sb = sum(Sb);
X = Y(:,delimiter(j):delimiter(j+1));
T = trace(Sb/Sw);
end

