function [coeff,linear_db,nonlinear_db] = mopso_fun(N_clu,n,Pvec,include_x,include_xp,include_crossterms,equalterms)
    %input: N class samples from the patient
    %output: after non linear trans coefficients
    linear_data  = DimReduc( N_clu,n,0);
    linear_db = linear_data;
 
    %linear data contains the patient's N data and V,S,F from DS1
    % data.N data.V. data.S data.F
    [data, nE] = nonlinear_trans(linear_data,n,Pvec,include_x,include_xp,include_crossterms,equalterms);

    %% Problem Definition
    CostFunction = @(x) CostFun(x,data); 
    nVar = nE;      % Number of Unknown (Decision) Variables
    VarMin = -1;   % Lower Bound of Decision Variables
    VarMax =  1;   % Upper Bound of Decision Variables

    VarSize=[1 nVar];   % Size of Decision Variables Matrix
    %% MOPSO Parameters

    MaxIt=150;           % Maximum Number of Iterations

    nPop=100;            % Population Size

    nRep=200;            % Repository Size

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
        %figure(1);
        %PlotCosts(pop,rep);
        %pause(0.01);

        % Show Iteration Information
        disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);

        % Damping Inertia Weight
        w=w*wdamp;

    end
    
    %% Results

    rep_costs=vertcat(rep.Cost);
    [~,I] = min(rep_costs(:,1) +rep_costs(:,2));
    coeff = rep(I).Position;
    nonlinear_db.N = data.N .*repmat(coeff,size(data.N,1),1);
    nonlinear_db.V = data.V .*repmat(coeff,size(data.V,1),1);
    nonlinear_db.S = data.S .*repmat(coeff,size(data.S,1),1);
    nonlinear_db.F = data.F .*repmat(coeff,size(data.F,1),1);
    labels = [2.*ones(size(nonlinear_db.V,1),1);...
        3.*ones(size(nonlinear_db.S,1),1); 4.*ones(size(nonlinear_db.F,1),1)];
    nonlinear_db.ab = [nonlinear_db.V; nonlinear_db.S; nonlinear_db.F];
    nonlinear_db.ab = [nonlinear_db.ab labels];
end
