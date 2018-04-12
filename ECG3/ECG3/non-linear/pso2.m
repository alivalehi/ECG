
clc;
clear;
close all;
%load('non_linear_1234.mat')
%load('cross_term_p3.mat')
%load('normlized_20.mat')
%load('matlab.mat', 'y')
%load('nMITDS1_4.mat')
load('\\EGRSHARES\Homes\NAU\jc3464\Documents\MATLAB\ECG2\MITDB(non-normalized)\100NF.mat')
linear_data  = DimReduc( Features(find(Features(:,29)==1),1:28),1);


%x_linear = nMITDS1_4(:,1:4); y = nMITDS1_4(:,5);clear nMITDS1_4;
n=4;Pvec=[1,2,3];	include_x = 0; include_xp = 0;include_crossterms = 1;equalterms = 0;
[ME, nE] = makePCoeffs(n,Pvec, include_x, include_xp, include_crossterms, equalterms);

if 0
XX = poly_tranform(x_linear,ME);
%inital = rand(1,12);
XX = zscore(XX);
%XX = tanh((XX-repmat(mean(XX),15357,1))./repmat(std(XX),15357,1));
data.N = XX(find(y==1),:);
data.V = XX(find(y==2),:);
data.S = XX(find(y==3),:);
data.F = XX(find(y==4),:);
end

data.N = poly_tranform(linear_data.N,ME);
data.V = poly_tranform(linear_data.V,ME);
data.S = poly_tranform(linear_data.S,ME);
data.F = poly_tranform(linear_data.F,ME);
%% Problem Definiton
data.all = zscore([data.N; data.V; data.S; data.F]);
data.N = data.all(1:size(data.N,1),:);
data.V = data.all(1+size(data.N,1):size(data.N,1)+size(data.V,1),:);
data.S = data.all(1+size(data.N,1)+size(data.V,1):size(data.N,1)+size(data.V,1)+size(data.S,1),:);
data.F = data.all(1+size(data.N,1)+size(data.V,1)+size(data.S,1):end,:);
problem.CostFunction = @(x) CostFun(x,data);  % Cost Function
problem.CostTerm = 5;
problem.nVar = nE;      % Number of Unknown (Decision) Variables
problem.VarMin = -10;   % Lower Bound of Decision Variables
problem.VarMax =  10;   % Upper Bound of Decision Variables

%% Parameters of PSO

% Constriction Coefficients
kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

params.MaxIt = 1000;        % Maximum Number of Iterations
params.nPop = 100;           % Population Size (Swarm Size)
params.w = chi;             % Intertia Coefficient
params.wdamp = 1;           % Damping Ratio of Inertia Coefficient
params.c1 = chi*phi1;       % Personal Acceleration Coefficient
params.c2 = chi*phi2;       % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

%% Calling PSO

out = PSO(problem, params);

BestSol = out.BestSol;
BestCosts = out.BestCosts;

%% Results

figure;
subplot(5,1,1)
semilogy(BestCosts(:,1), 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
subplot(5,1,2)
semilogy(BestCosts(:,2), 'LineWidth', 2);
xlabel('Iteration');
ylabel('term 1');
grid on;
subplot(5,1,3)
semilogy(BestCosts(:,3), 'LineWidth', 2);
xlabel('Iteration');
ylabel('term 2');
grid on;
subplot(5,1,4)
semilogy(BestCosts(:,4), 'LineWidth', 2);
xlabel('Iteration');
ylabel('term 3');
grid on;
subplot(5,1,5)
semilogy(BestCosts(:,5), 'LineWidth', 2);
xlabel('Iteration');
ylabel('term 4');
grid on;

test_data = [data.N;data.V;data.S;data.F];
tt = test_data(:,1:nE).*repmat(out.BestSol.Position,size(test_data,1),1);
clear test_data
test_data = [tt [ones(size(data.N,1),1);2.*ones(size(data.V,1),1);3.*ones(size(data.S,1),1);4.*ones(size(data.F,1),1)]];