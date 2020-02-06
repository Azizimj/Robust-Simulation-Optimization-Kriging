function out = PSO(problem, params)
global x_mat %simulation points n*d
global w_mean
global e_mat
global n
global x_pred
    %% Problem Definiton

    CostFunction = problem.CostFunction;  % Cost Function

    nvar = problem.nVar;        % Number of Unknown (Decision) Variables
    mvar = problem.mVar;        % Number of noise Variables
    VarSize = [1 nvar];         % Matrix Size of Decision Variables

    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables
    e_high=problem.NoiseMax;
    e_low=problem.NoiseMin;
    

    %% Parameters of PSO

    PsoMaxIt = params.PsoMaxIt;   % Maximum Number of Iterations for PSO
    AlgoMaxIt= params.AlgoMaxIt;   % Maximum Number of Iterations for algorithm
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
    empty_particle.noise=[];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.noise = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost = inf;
    GlobalBest.Position = [];
    GlobalBest.noise = [];
    GlobalBest.Position=[];

   %% Alghorithm
    qual=inf; % solution quality
    thersh=.001; 
    counter=0;
    BestCosts = zeros(PsoMaxIt*AlgoMaxIt, 1);
    lhsused=0;
    while (qual >= thersh) && (counter<AlgoMaxIt)   
    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Initialize Velocity
        particle(i).Velocity = unifrnd(zeros(1,nvar),0.2*(VarMax-VarMin),VarSize);
        
        % Evaluation
%         e=sym('e',[1 mvar]);
%         opts = optimset('Display','iter');
%         [ee,fval]= fmincon(matlabFunction(-1*CostFunction(particle(i).Position,e)),ones(1,nvar),[],[],[],[],low_noise,high_noise);
%         particle(i).Cost = -1*fval;
%         s=solve(diff(CostFunction(particle(i).Position,e))==0,e);
%         particle(i).Cost = CostFunction(particle(i).Position,s);
%       
%         particle(i).Cost = CostFunction(particle(i).Position,ones(1,mvar));
%         particle(i).noise= ones(1,mvar);
        
        x_pred=particle(i).Position;
        [ee,fval]= fmincon(@kriging,ones(1,mvar),[],[],[],[],e_low,e_high);
        eee(i,:)=ee;
        particle(i).Cost = -1*fval;
        particle(i).noise= ee;
        
        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.noise=particle(i).noise;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest.Position = particle(i).Best.Position;
            GlobalBest.noise= particle(i).Best.noise;
            GlobalBest.Cost = particle(i).Best.Cost;
        end

    end

    % Array to Hold Best Cost Value on Each Iteration
   


    %% Main Loop of PSO

    for it=1:PsoMaxIt

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
%             particle(i).Cost = CostFunction(particle(i).Position);
%             e=sym('e',[1 mvar]);
%             [ee,fval2]=fmincon(matlabFunction(-1*CostFunction(particle(i).Position,e)),zeros(1,nvar),[],[],[],[],low_noise,high_noise);             
%             particle(i).Cost = -1*fval2;
%             particle(i).Cost = CostFunction(particle(i).Position,ones(1,mvar));
%             particle(i).noise= ones(1,mvar);
            
        x_pred=particle(i).Position;
        [ee,fval]= fmincon(@kriging,ones(1,mvar),[],[],[],[],e_low,e_high);
        eee1(i+(it-1)*nPop,:)=ee;
        particle(i).Cost = -1*fval;
        particle(i).noise= ee;
        
            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.noise = particle(i).noise;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end            

            end

        end

        % Store the Best Cost Value
        BestCosts(counter*PsoMaxIt+it) = GlobalBest.Cost;

        % Display Iteration Information
        if ShowIterInfo
            disp(['AlgoIteration ' num2str(counter+1) ' PSoIteration ' num2str(it)...
                ': Best Cost = ' num2str(BestCosts(counter*PsoMaxIt+it))]);
        end

        % Damping Inertia Coefficient
        w = w * wdamp;

    end
    sim=simulation(GlobalBest.Position,GlobalBest.noise);
    qual= 2*abs(sim-GlobalBest.Cost)/(abs(sim)+abs(GlobalBest.Cost)); 
         if (qual>=thersh)
            lhsused=lhsused+1;
            x_r=(VarMax-VarMin)/2; e_r=(e_high-e_low)/2;
            damp=.95;
            L=LHS(max(GlobalBest.Position-(x_r)*damp,VarMin),min(GlobalBest.Position+(x_r)*damp,VarMax),...
                 max(GlobalBest.noise-(e_r)*damp,e_low),min(GlobalBest.noise+(e_r)*damp,e_high),n);
            VarMin=max(GlobalBest.Position-(x_r)*damp,VarMin);
            VarMax=min(GlobalBest.Position+(x_r)*damp,VarMax);
            e_low=max(GlobalBest.noise-(e_r)*damp,e_low);
            e_high=min(GlobalBest.noise+(e_r)*damp,e_high);
            MaxVelocity = 0.2*(VarMax-VarMin);
            MinVelocity = -MaxVelocity;
         
            x_mat=L.x;
            e_mat=L.e;
        
            for i=1:size(x_mat,1)
                sim=simulation(x_mat(i,:),e_mat(i,:));
                w_mean(i)=sim;
            end
        end
        counter=counter+1;
    end
    
    out.pop = particle;
    out.BestSol = GlobalBest;
    out.BestCosts = BestCosts;
    out.lhs=lhsused;
end