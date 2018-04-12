function out = PSO(problem, params)

    %% Problem Definiton

    CostFunction = problem.CostFunction;  % Cost Function
    
    cost_term = problem.CostTerm;
    
    data = problem.data;

    nVar = problem.nVar;        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];         % Matrix Size of Decision Variables

    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables


    %% Parameters of PSO

    MaxIt = params.MaxIt;   % Maximum Number of Iterations

    nPop = params.nPop;     % Population Size (Swarm Size)

    w = params.w;           % Intertia Coefficient
    wdamp = params.wdamp;   % Damping Ratio of Inertia Coefficient
    c1 = params.c1;         % Personal Acceleration Coefficient
    c2 = params.c2;         % Social Acceleration Coefficient

    % The Flag for Showing Iteration Information
    ShowIterInfo = params.ShowIterInfo;    

    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    
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
    GlobalBest.Cost.Term = NaN(1,5);

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        [particle(i).Cost.Sum,particle(i).Cost.Term] = CostFunction(particle(i).Position,data);

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost.Sum < GlobalBest.Cost.Sum
%         score = sum([particle(i).Best.Cost.Term(1) < GlobalBest.Cost.Term(1),...
%             particle(i).Best.Cost.Term(2) < GlobalBest.Cost.Term(2),...
%             particle(i).Best.Cost.Term(3) < GlobalBest.Cost.Term(3)]);
        %if particle(i).Best.Cost.Term(1) < GlobalBest.Cost.Term(1) && particle(i).Best.Cost.Term(2) < GlobalBest.Cost.Term(2)...
        %        &&particle(i).Best.Cost.Term(3) < GlobalBest.Cost.Term(3)
        %if score >= 2
            GlobalBest = particle(i).Best;
            disp('found global');
        end

    end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt, cost_term+1);


    %% Main Lo1 of PSO

    for it=1:MaxIt

        for i=1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            % Evaluation
            [particle(i).Cost.Sum,particle(i).Cost.Term] = CostFunction(particle(i).Position,data);

            % Update Personal Best
            if particle(i).Cost.Sum < particle(i).Best.Cost.Sum
%             if particle(i).Cost.Term(1) < particle(i).Best.Cost.Term(1) && particle(i).Cost.Term(2) < particle(i).Best.Cost.Term(2)...
%                     &&particle(i).Best.Cost.Term(3) < GlobalBest.Cost.Term(3)
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                 if particle(i).Best.Cost.Sum < GlobalBest.Cost.Sum
%                 if particle(i).Best.Cost.Term(1) < GlobalBest.Cost.Term(1) && particle(i).Best.Cost.Term(2) < GlobalBest.Cost.Term(2)...
%                         &&particle(i).Best.Cost.Term(3) < GlobalBest.Cost.Term(3)
                    GlobalBest = particle(i).Best;
                end            

            end

        end

        % Store the Best Cost Value
        BestCosts(it,1) = GlobalBest.Cost.Sum;
        BestCosts(it,2:end) = GlobalBest.Cost.Term;
        

        % Display Iteration Information
        if ShowIterInfo
            disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it,1))]);
        end

        % Damping Inertia Coefficient
        w = w * wdamp;

    end
    
    out.pop = particle;
    out.BestSol = GlobalBest;
    out.BestCosts = BestCosts;
    
end