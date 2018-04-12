%%%%%% dummy data test
clear
addpath('non-linear/MOPSO')
% initialize non-linear parameters
n_lin=2;Pvec=[1,2];	include_x = 1; include_xp = 1;include_crossterms = 1;equalterms = 1;

mu1 = [1,1];
sigma1 = [1,0.5];
rng default
r1 = mvnrnd(mu1,sigma1,1000);
lin_data.N = r1;

mu2 = [2,5];
sigma2 = [0.5,1];
rng default
r2 = mvnrnd(mu2,sigma2,300);
lin_data.S = r2;

mu3 = [6,8];
sigma3 = [1,0.5];
rng default
r3 = mvnrnd(mu3,sigma3,500);
lin_data.V = r3;

mu4 = [8,3];
sigma4 = [1,1];
rng default
r4 = mvnrnd(mu4,sigma4,100);
lin_data.F = r4;

figure(1);scatter(r1(:,1),r1(:,2),'bd');
hold on;scatter(r2(:,1),r2(:,2),'r+');
hold on;scatter(r3(:,1),r3(:,2),'co');
hold on;scatter(r4(:,1),r4(:,2),'k*');
title('dummy 2-D data','FontSize',14)
hold off;
x= 5,y=5,r=2;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
lin_data = [xunit',yunit'];
%figure;
%plot(xunit, yunit);
 [data, nE] = nonlinear_trans(lin_data,n_lin,Pvec,include_x,include_xp,include_crossterms,equalterms);
 ali=1;
%  r5=data.N; 
%  r6= data.S;
%  r7=data.V;
%  r8=data.F;
 i=3;
 j=3;
%figure(1);scatter(data(:,i),data(:,j)+data(:,j+1),'bd');
%  figure(1);scatter(r5(:,i),r5(:,j),'bd');
% hold on;scatter(r6(:,i),r6(:,j),'r+');
% hold on;scatter(r7(:,i),r7(:,j),'co');
% hold on;scatter(r8(:,i),r8(:,j),'k*');
% title('dummy 2-D data','FontSize',14)
% hold off;
% %[data, nE] = transform_x(lin_data);
% 
% %% Problem Definition
% CostFunction = @(x) CostFun(x,data); 
% nVar = nE;      % Number of Unknown (Decision) Variables
% VarMin = -1;   % Lower Bound of Decision Variables
% VarMax =  1;   % Upper Bound of Decision Variables
% 
% VarSize=[1 nVar];   % Size of Decision Variables Matrix
% %% MOPSO Parameters
% 
% MaxIt=150;           % Maximum Number of Iterations
% 
% nPop=100;            % Population Size
% 
% nRep=200;            % Repository Size
% 
% w=0.5;              % Inertia Weight
% wdamp=0.99;         % Intertia Weight Damping Rate
% c1=1;               % Personal Learning Coefficient
% c2=2;               % Global Learning Coefficient
% 
% nGrid=7;            % Number of Grids per Dimension
% alpha=0.1;          % Inflation Rate
% 
% beta=2;             % Leader Selection Pressure
% gamma=2;            % Deletion Selection Pressure
% 
% mu=0.1;             % Mutation Rate
% 
% %% Initialization
% 
% empty_particle.Position=[];
% empty_particle.Velocity=[];
% empty_particle.Cost=[];
% empty_particle.Best.Position=[];
% empty_particle.Best.Cost=[];
% empty_particle.IsDominated=[];
% empty_particle.GridIndex=[];
% empty_particle.GridSubIndex=[];
% 
% pop=repmat(empty_particle,nPop,1);
% 
% for i=1:nPop
% 
%     pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
% 
%     pop(i).Velocity=zeros(VarSize);
% 
%     pop(i).Cost=CostFunction(pop(i).Position);
% 
% 
%     % Update Personal Best
%     pop(i).Best.Position=pop(i).Position;
%     pop(i).Best.Cost=pop(i).Cost;
% 
% end
% 
% % Determine Domination
% pop=DetermineDomination(pop);
% 
% rep=pop(~[pop.IsDominated]);
% 
% Grid=CreateGrid(rep,nGrid,alpha);
% 
% for i=1:numel(rep)
%     rep(i)=FindGridIndex(rep(i),Grid);
% end
% 
% 
% %% MOPSO Main Loop
% 
% for it=1:MaxIt
% 
%     for i=1:nPop
% 
%         leader=SelectLeader(rep,beta);
% 
%         pop(i).Velocity = w*pop(i).Velocity ...
%             +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
%             +c2*rand(VarSize).*(leader.Position-pop(i).Position);
% 
%         pop(i).Position = pop(i).Position + pop(i).Velocity;
% 
%         pop(i).Position = max(pop(i).Position, VarMin);
%         pop(i).Position = min(pop(i).Position, VarMax);
% 
%         pop(i).Cost = CostFunction(pop(i).Position);
% 
%         % Apply Mutation
%         pm=(1-(it-1)/(MaxIt-1))^(1/mu);
%         if rand<pm
%             NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
%             NewSol.Cost=CostFunction(NewSol.Position);
%             if Dominates(NewSol,pop(i))
%                 pop(i).Position=NewSol.Position;
%                 pop(i).Cost=NewSol.Cost;
% 
%             elseif Dominates(pop(i),NewSol)
%                 % Do Nothing
% 
%             else
%                 if rand<0.5
%                     pop(i).Position=NewSol.Position;
%                     pop(i).Cost=NewSol.Cost;
%                 end
%             end
%         end
% 
%         if Dominates(pop(i),pop(i).Best)
%             pop(i).Best.Position=pop(i).Position;
%             pop(i).Best.Cost=pop(i).Cost;
% 
%         elseif Dominates(pop(i).Best,pop(i))
%             % Do Nothing
% 
%         else
%             if rand<0.5
%                 pop(i).Best.Position=pop(i).Position;
%                 pop(i).Best.Cost=pop(i).Cost;
%             end
%         end
% 
%     end
% 
%     % Add Non-Dominated Particles to REPOSITORY
%     rep=[rep
%          pop(~[pop.IsDominated])]; %#ok
% 
%     % Determine Domination of New Resository Members
%     rep=DetermineDomination(rep);
% 
%     % Keep only Non-Dminated Memebrs in the Repository
%     rep=rep(~[rep.IsDominated]);
% 
%     % Update Grid
%     Grid=CreateGrid(rep,nGrid,alpha);
% 
%     % Update Grid Indices
%     for i=1:numel(rep)
%         rep(i)=FindGridIndex(rep(i),Grid);
%     end
% 
%     % Check if Repository is Full
%     if numel(rep)>nRep
% 
%         Extra=numel(rep)-nRep;
%         for e=1:Extra
%             rep=DeleteOneRepMemebr(rep,gamma);
%         end
% 
%     end
% 
%     % Plot Costs
%     figure(5);
%     PlotCosts(pop,rep);
%     pause(0.01);
% 
%     % Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
% 
%     % Damping Inertia Weight
%     w=w*wdamp;
% 
% end
% 
% %% Results
% theta = 0.1;
% rep_costs=vertcat(rep.Cost);
% [~,I] = min(theta*rep_costs(:,1) + (1-theta)*rep_costs(:,2));
% coeff = rep(I).Position;
% nonlinear_db.N = data.N .*repmat(coeff,size(data.N,1),1);
% nonlinear_db.V = data.V .*repmat(coeff,size(data.V,1),1);
% nonlinear_db.S = data.S .*repmat(coeff,size(data.S,1),1);
% nonlinear_db.F = data.F .*repmat(coeff,size(data.F,1),1);
% labels = [2.*ones(size(nonlinear_db.V,1),1);...
%     3.*ones(size(nonlinear_db.S,1),1); 4.*ones(size(nonlinear_db.F,1),1)];
% nonlinear_db.ab = [nonlinear_db.V; nonlinear_db.S; nonlinear_db.F];
% nonlinear_db.ab = [nonlinear_db.ab labels];
% nonlinear_db.all = [nonlinear_db.N ones(size(nonlinear_db.N,1),1);nonlinear_db.ab];
% all = nonlinear_db.all;
% 
% [coeff_pca,score,latent,tsquared,explained,mu] = pca(all(:,1:end-1));
% visual = all(:,1:end-1)*coeff_pca(:,1:2);
% visual = [visual all(:,end)];
% 
% figure(2);scatter(visual(all(:,end)==1,1),visual(all(:,end)==1,2),'bd');
% hold on;scatter(visual(all(:,end)==3,1),visual(all(:,end)==3,2),'r+');
% hold on;scatter(visual(all(:,end)==2,1),visual(all(:,end)==2,2),'co');
% hold on;scatter(visual(all(:,end)==4,1),visual(all(:,end)==4,2),'*');
% hold off;
% 
% figure(3)
% %gplotmatrix(all(:,1:end-1),[],all(:,end),['b','r' 'c' 'k'],['d' '+' 'o' '*'],[],false);
% scatter3(all(all(:,4)==1,1),all(all(:,4)==1,2),all(all(:,4)==1,3),'bd');
% hold on;scatter3(all(all(:,4)==2,1),all(all(:,4)==2,2),all(all(:,4)==2,3),'r+');
% hold on;scatter3(all(all(:,4)==3,1),all(all(:,4)==3,2),all(all(:,4)==3,3),'co');
% hold on;scatter3(all(all(:,4)==4,1),all(all(:,4)==4,2),all(all(:,4)==4,3),'k*');
% title('transformed data','FontSize',14)
% hold off;
% 
