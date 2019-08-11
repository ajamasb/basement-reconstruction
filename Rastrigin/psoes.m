function [varargout] = psoes(CostFunction, VarSize, VarMin, VarMax, varargin)
% To be written ...

narginchk(3, inf);


if (numel(VarMin) ~= numel(VarMax) || any(VarMin>= VarMax))
    error('pso:ErrorVar',['Error.\nThe bounds of the search space are not consistent\n',...
        'This error may be caused by either unequal number of upper and lower bounds\n',...
        'or in any case where an element of Varmin is greater than any elements of VarMax'])

    
end

%% input Parser
p = inputParser;

%General PSO information

defaultMaxIt=100;      
defaultnPop=25;  
defaulDisplay = true;

% PSO Parameters (PSO with inertia weight)

defaultW=1;            
defaultWdamp=0.99;     
defaultC1=1.5;         
defaultC2=2.0;         

% PSO Paramaters (PSO with constriction factor)

defaultConstrictionFactor = false;
defaultPhi1=2.05;
defaultPhi2=2.05;

% PSOES parameters (mutation)
defaulMutation = true(1);
defaultc_plus=1.75;
defaultc_minus=0.97;
defaultmutLimit=5;



%%

% General PSO
addRequired(p, 'CostFunction', @(x) (isa(x,'function_handle')))
addRequired(p, 'VarSize', @(x) (validateattributes(x,{'numeric'},{'numel',2})))
addRequired(p, 'VarMin', @(x) (validateattributes(x,{'numeric'},{'2d','nonempty'})))
addRequired(p, 'VarMax', @(x) (validateattributes(x,{'numeric'},{'2d','nonempty'})))
addParameter(p, 'nPop', defaultnPop, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'MaxIt', defaultMaxIt, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'Display',defaulDisplay,@(x) (validateattributes(x, {'logical'}, {'nonempty'})))

% PSO with Inertia Weight
addParameter(p, 'w', defaultW, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'c1', defaultC1, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'c2', defaultC2, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))

% PSO with Constriction Factor
addParameter(p, 'ConstrictionFactor', defaultConstrictionFactor,@(x) (validateattributes(x, {'logical'}, {'nonempty'})))
addParameter(p, 'phi1', defaultPhi1, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'phi2', defaultPhi2, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))

% PSOES: PSO with mutation
addParameter(p, 'mutation',defaulMutation,@(x) (validateattributes(x, {'logical'}, {'nonempty'})))
addParameter(p, 'c_minus', defaultc_minus, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'c_plus', defaultc_plus, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))
addParameter(p, 'MuteLimit', defaultmutLimit, @(x) (validateattributes(x, {'numeric'}, {'scalar'})))



parse (p, CostFunction, VarSize, VarMin, VarMax, varargin{:});


%% Mutation Constants
step=1;
cmut=0;
ismutated=false;



%======================================================================================
if p.Results.ConstrictionFactor
    
    phi= p.Results.phi1+ p.Results.phi2;
    chi=2/(phi-2+sqrt(phi^2-4*phi));
    w=chi;          % Inertia Weight
    wdamp=1;        % Inertia Weight Damping Ratio
    c1=chi* p.Results.phi1;    % Personal Learning Coefficient
    c2=chi* p.Results.phi2;    % Global Learning Coefficient
    
else
    
    w =  p.Results.w;
    c1 =  p.Results.c1;
    c2 =  p.Results.c2;
    wdamp = 0.99;
    
end

%% Velocity Limits 


VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Creating the swarm (here named particle)

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,p.Results.nPop,1);

GlobalBest.Cost=inf;
% [GlobalBest,particle] = ScatterParticles( particle, CostFunction, VarMin, VarMax, VarSize,p.Results.nPop);
for i=1:p.Results.nPop
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position); 
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end


%% main pso loop
BestCost = zeros(p.Results.MaxIt,1);
for it=1:p.Results.MaxIt
    
    for i=1:p.Results.nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        
%       % Mutation 
        if p.Results.mutation
        
         
            if rand<1;
                PosMutated.x = mutate(particle(i).Position,VarSize(2),step,VarMin,VarMax);
                
                if (cmut) <p.Results.MuteLimit;
                    step = step * p.Results.c_minus;
                else
                    step = step * p.Results.c_plus;
                    cmut=0;
                end
                
                
                ismutated=true;
                
            end  
               if (ismutated)
                   PosMutated.cost = CostFunction(PosMutated.x);
                   
                   
                   if (PosMutated.cost<particle(i).Best.Cost)
                       particle(i).Position = PosMutated.x;
                       particle(i).Cost = PosMutated.cost; 
                       cmut=cmut+1;
                       
                       ismutated = false;
                                                 
                   end
               
                                 
                   
               end 
        end
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
  
    if p.Results.Display ;disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) ' Step: ' num2str(step)]);end
    
    w=w*wdamp;
    

  
end

varargout{1} = GlobalBest;
varargout{2} = BestCost;
varargout{3} = p.Results;
end







