% non-linear tranformation test 4-12-2017
clear
load('nMITDS1_4.mat')
x_linear = nMITDS1_4(:,1:4); y = nMITDS1_4(:,5);clear nMITDS1_4;
n=4;Pvec=[1,2,3];	include_x = 0; include_xp = 1;include_crossterms = 1;equalterms = 1;
[ME, nE] = makePCoeffs(n,Pvec, include_x, include_xp, include_crossterms, equalterms);
XX = poly_tranform(x_linear,ME);
%inital = rand(1,12);
XX = zscore(XX);
%XX = tanh((XX-repmat(mean(XX),15357,1))./repmat(std(XX),15357,1));
data.N = XX(find(y==1),:);
data.V = XX(find(y==2),:);
data.S = XX(find(y==3),:);
data.F = XX(find(y==4),:);
