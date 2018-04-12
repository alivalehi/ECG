% generate toydata
clear
mu1 = [0 1 2 3];
sigma1 = [1 0 0 0;0 2 0 0; 0 0 3 0;0 0 0 1];
%sigma1 = [0.01 0 0;0 0.02 0; 0 0 .03];
rng default
c1 = mvnrnd(mu1,sigma1,300);

mu2 = [2 -2 0 -1];
%sigma2 = [1 .5 0;.5 2 .5;0 .5 3];
sigma2 = [1 0 0 0;0 2 0 0; 0 0 3 0;0 0 0 1];%[.01 .005 0 0;.005 .02 .005 0;0 .005 .03 0; 0 0 0 0.5];
rng default
c2 = mvnrnd(mu2,sigma2,500);

mu3 = [1 .5 1 -2];
sigma3 = [1 0 0 0;0 2 0 0; 0 0 3 0;0 0 0 1];%[1 .5 0 0;.5 2 .5 0;0 .5 3 0; 0 0 0 0.8];
rng default
c3 = mvnrnd(mu3,sigma3,700);

mu4 = [2 -1 0 0];
sigma4 = [1 0 0 0;0 2 0 0; 0 0 3 0;0 0 0 1];%[1 .5 0 0.1;.5 2 .5 0;0 .5 3 0; 0.1 0 0 .9];
rng default
c4 = mvnrnd(mu4,sigma4,200);

data.all = zscore([c2;c3;c1;c4]);
data.N = data.all(1:size(c2,1),:);
data.V = data.all(1+size(c2,1):size(c2,1)+size(c3,1),:);
data.S = data.all(1+size(c2,1)+size(c3,1):size(c2,1)+size(c3,1)+size(c1,1),:);
data.F = data.all(1+size(c2,1)+size(c3,1)+size(c1,1):end,:);

save('toydata.mat','data');