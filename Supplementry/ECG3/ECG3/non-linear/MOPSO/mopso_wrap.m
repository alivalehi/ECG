clc;
clear;
close all;

test = 0;
%% Loading Data
switch test
    case 0
        load('\\EGRSHARES\Homes\NAU\jc3464\Documents\MATLAB\ECG2\MITDB(non-normalized)\105NF.mat')
        linear_data  = DimReduc( Features(find(Features(:,29)==1),1:28),1);

    case 1
        load('toydata.mat');
        linear_data = data;
        clear data;
end

%x_linear = nMITDS1_4(:,1:4); y = nMITDS1_4(:,5);clear nMITDS1_4;
n=4;Pvec=[1,2,3];	include_x = 1; include_xp = 1;include_crossterms = 1;equalterms = 1;
[ME, nE] = makePCoeffs(n,Pvec, include_x, include_xp, include_crossterms, equalterms);

data.N = poly_tranform(linear_data.N,ME);
data.V = poly_tranform(linear_data.V,ME);
data.S = poly_tranform(linear_data.S,ME);
data.F = poly_tranform(linear_data.F,ME);

data.all = zscore([data.N; data.V; data.S; data.F]);
data.N = data.all(1:size(data.N,1),:);
data.V = data.all(1+size(data.N,1):size(data.N,1)+size(data.V,1),:);
data.S = data.all(1+size(data.N,1)+size(data.V,1):size(data.N,1)+size(data.V,1)+size(data.S,1),:);
data.F = data.all(1+size(data.N,1)+size(data.V,1)+size(data.S,1):end,:);
%% Problem Definition
CostFunction = @(x) CostFun(x,data); 
nVar = nE;      % Number of Unknown (Decision) Variables
VarMin = -10;   % Lower Bound of Decision Variables
VarMax =  10;   % Upper Bound of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

%% MOPSO Parameters

MaxIt=200;           % Maximum Number of Iterations

nPop=200;            % Population Size

nRep=100;            % Repository Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

nGrid=7;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate

%% Initialization

empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];

pop=repmat(empty_particle,nPop,1);

for i=1:nPop
    
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Velocity=zeros(VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    
    % Update Personal Best
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
    
end

% Determine Domination
pop=DetermineDomination(pop);

rep=pop(~[pop.IsDominated]);

Grid=CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end


%% MOPSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        leader=SelectLeader(rep,beta);
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Position = max(pop(i).Position, VarMin);
        pop(i).Position = min(pop(i).Position, VarMax);
        
        pop(i).Cost = CostFunction(pop(i).Position);
        
        % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
            NewSol.Cost=CostFunction(NewSol.Position);
            if Dominates(NewSol,pop(i))
                pop(i).Position=NewSol.Position;
                pop(i).Cost=NewSol.Cost;

            elseif Dominates(pop(i),NewSol)
                % Do Nothing

            else
                if rand<0.5
                    pop(i).Position=NewSol.Position;
                    pop(i).Cost=NewSol.Cost;
                end
            end
        end
        
        if Dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            
        elseif Dominates(pop(i).Best,pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
    figure(1);
    PlotCosts(pop,rep);
    pause(0.01);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w=w*wdamp;
    
end

%% Resluts

rep_costs=vertcat(rep.Cost);
thres = prctile(rep_costs,20);
cand = [];
% for i=1:length(rep_costs)
%     if rep_costs(i,1)<thres(1) && rep_costs(i,2)<thres(2)
%         cand(end+1) = i;
%     end
% end
[M,I] = min(rep_costs(:,1) +rep_costs(:,2));
%rep(I).Cost

test_data = [data.N;data.V;data.S;data.F];
tt = test_data(:,1:nE).*repmat(rep(I).Position,size(test_data,1),1);
clear test_data
test_data = [tt [ones(size(data.N,1),1);2.*ones(size(data.V,1),1);3.*ones(size(data.S,1),1);4.*ones(size(data.F,1),1)]];