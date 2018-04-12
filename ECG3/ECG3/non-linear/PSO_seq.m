function out = PSO_seq(problem, params)

    %% Problem Definiton

    CostFunction = problem.CostFunction;  % Cost Function
    
    cost_term = problem.CostTerm;       % The total number of cost terms that we want to control
    
    data = problem.data;

    nVar = problem.nVar;        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];
    % Matrix Size of Decision Variables

    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables
    
    op = problem.OptimizeOver;
    sl = problem.SelectOver;


    %% Parameters of PSO

    MaxIt = params.MaxIt;   % Maximum Number of Iterations

    nPop = params.nPop;     % Population Size (Swarm Size)
    
    
    w = params.w;           % Inertia Coefficient
    wdamp = params.wdamp;   % Damping Ratio of Inertia Coefficient
    c1 = params.c1;         % Personal Acceleration Coefficient
    c2 = params.c2;         % Social Acceleration Coefficient

    % The Flag for Showing Iteration Information
    ShowIterInfo = params.ShowIterInfo;    

    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    
    %% parameters for system
    % Mask of parameter, initialized as all zeros
    Mask = logical([zeros(1,nVar)]);
    if problem.IncludeOrigin
        ParticleSize = [1 4];
        for x = 1:4
            Mask(x) = true;
        end
        start = 4;
    else
        ParticleSize = [1 1];
        start = 1;
    end
           % we define one dimension particles to optimize one variable
    % For each round we optimize one variable according to cost term 1 and check the final value of

    
for v = start:nVar
    % we're optimizing parameter v now
    Mask(v) = true;
    Masked_data.N = data.N(:,Mask);
    Masked_data.V = data.V(:,Mask);
    Masked_data.S = data.S(:,Mask);
    Masked_data.F = data.F(:,Mask);
    %% Initialization

    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost.Sum = inf;
    GlobalBest.Cost.Term = inf(1,5);

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = zscore(unifrnd(VarMin, VarMax, ParticleSize));

        % Initialize Velocity
        % particle(i).Velocity = zeros(VarSize);
        particle(i).Velocity = zeros(ParticleSize);

        % Evaluation
        [particle(i).Cost.Sum,particle(i).Cost.Term] = CostFunction(particle(i).Position,Masked_data);

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost.Term(op) < GlobalBest.Cost.Term(op) && particle(i).Best.Cost.Term(2) < GlobalBest.Cost.Term(2)
            GlobalBest = particle(i).Best;
        end
    end


    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt, cost_term+1);


    %% Main Loop of PSO

    for it=1:MaxIt

        for i=1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(ParticleSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(ParticleSize).*(GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            % Evaluation
            [particle(i).Cost.Sum,particle(i).Cost.Term] = CostFunction(particle(i).Position,Masked_data);

            % Update Personal Best
            if particle(i).Cost.Term(op) < particle(i).Best.Cost.Term(op) && particle(i).Cost.Term(2) < particle(i).Best.Cost.Term(2)

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost.Term(op) < GlobalBest.Cost.Term(op) && particle(i).Best.Cost.Term(2) < GlobalBest.Cost.Term(2)
                    GlobalBest = particle(i).Best;
                end            

            end

        end

        % Store the Best Cost Value
        BestCosts(it,1) = GlobalBest.Cost.Sum;
        BestCosts(it,2:end) = GlobalBest.Cost.Term;
        

        % Damping Inertia Coefficient
        w = w * wdamp;

    end
    
    out(v-start+1).pop = particle;
    out(v-start+1).BestSol = GlobalBest;
    out(v-start+1).BestCosts = BestCosts;
    out(v-start+1).Mask = Mask;
    if v == start
        ParticleSize(2)  = ParticleSize(2)+1;
        DimOptCost = out(v-start+1).BestCosts(MaxIt,sl+1);
        continue
    else
        %if out(v-start+1-1).BestCosts(MaxIt,1)<= out(v-start+1).BestCosts(MaxIt,1)
        %if out(v-start+1-1).BestCosts(MaxIt,3)>=1
        if out(v-start+1).BestCosts(MaxIt,sl+1)> DimOptCost
            Mask(v)= false;
%         elseif out(v-start+1-1).BestCosts(MaxIt,3)<= out(v-start+1).BestCosts(MaxIt,3)
%              Mask(v)= false;   
        else
            ParticleSize(2)  = ParticleSize(2)+1;
            DimOptCost = out(v-start+1).BestCosts(MaxIt,sl+1);
        end
        
        % Display Iteration Information
        if ShowIterInfo
            disp(['Iteration ' num2str(it) ': Term 1 = ' num2str(DimOptCost)]);
        end
    end
    
end
end