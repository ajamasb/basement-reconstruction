clc;clear;close all

%% Define the search-space
addpath('.\Src\')
load SYNTHETIC.mat

[VarMin, VarMax, nVar, VarSize, data] = SearchSpace(SYNTHETIC);



%% defining the cost function
CostFunction = @(x)(Cost2Dgzz(x,data));



%% perform the inversion - tune the PSOES parameters in optional pair values:

[sol,sol.ConvHis] = psoes(CostFunction,VarSize,VarMin,VarMax,'MaxIt',500,'c2',1.5,'nPop',30,'c_minus',0.95);

%% Plotting the results
figure

plot(data.obs.xobs,-sol.Position,'LineWidth',3)
hold on
plot(SYNTHETIC.Obs.xobs,SYNTHETIC.TrueModel,'LineWidth',3)
plot(SYNTHETIC.Obs.xobs,-VarMax,'m','LineWidth',3)
plot(SYNTHETIC.Obs.xobs,-VarMin,'m','LineWidth',3)
grid minor
xlabel('Distance')
ylabel('Depth (m)')
title('Rersults of the Inversion')
legend('RECONSTRUCTED','TRUE','UPPER BOUND OF SS','LOWER BOUND OF SS')


figure('color', [ 1 1 1])
plot(1:500,sol.ConvHis,'k', 'LineWidth',3)
axis tight
grid minor
title('CONVERGENCE HISTORY')
xlabel('Iteration')
ylabel('Cost')
set(gca,'fontweight','bold')
