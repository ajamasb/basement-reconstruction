%% producing figure 3 


clc;clear;close all

[VarMin, VarMax, nVar, VarSize] = SearchSpace();
A=10; CostFunction=@(x) (2*A+ (x(1).^2-A*cos(2*pi*x(1)))+(x(2).^2-A*cos(2*pi*x(2)))); D=2;



MaxRun = 100;

sol= cell(MaxRun+1,3);
sol{1} = 'PSOES';
sol{1,2} = 'PSO with inertia weight';
sol{1,3} = 'PSO with constriction factor';

Cost = zeros(MaxRun,3);

for run = 1: MaxRun
    
    
    
    
    sol{run+1,1} = psoes(CostFunction,VarSize,VarMin,VarMax,'Display',false);
    sol{run+1,2} = psoes(CostFunction,VarSize,VarMin,VarMax,'Mutation',false,'Display',false);
    sol{run+1,3} = psoes(CostFunction,VarSize,VarMin,VarMax,'Mutation',false, 'ConstrictionFactor',true,'Display',false);
    fprintf('Run: %d cost:   \n psoes: %g - pso inertia : %g - pso Constriction: %g\n',...
        run,sol{run+1,1}.Cost,sol{run+1,2}.Cost,sol{run+1,3}.Cost)
    Cost(run,:) = [sol{run+1,1}.Cost,sol{run+1,2}.Cost,sol{run+1,3}.Cost];
end


%%

figure('color',[1 1 1],'Position',[307.4000 430.6000 951.2000 322.4000])

subplot(1,3,1)
histogram(log10(Cost(:,1)),'facecolor','r')
axis tight
title('PSOES')
xlabel('log(Global Best Cost)')
set(gca,'fontweight','bold')

subplot(1,3,2)
histogram(log10(Cost(:,2)),'facecolor','r')
axis tight
title('PSO Inertia')
xlabel('log(Global Best Cost)')
set(gca,'fontweight','bold')

subplot(1,3,3)
histogram(log10(Cost(:,3)),'facecolor','r')
axis tight
title('PSO Constriction')
xlabel('log(Global Best Cost)')
set(gca,'fontweight','bold')


%%
format long
min(Cost)