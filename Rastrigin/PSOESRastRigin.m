clc;clear;close all



%% Defining the Cost Function

A=10; CostFunction=@(x) (2*A+ (x(1).^2-A*cos(2*pi*x(1)))+(x(2).^2-A*cos(2*pi*x(2)))); D=2;

% Again for plotting:
A=10; Rast=@(x,y) (2*A+ (x.^2-A*cos(2*pi*x))+(y.^2-A*cos(2*pi*y))); D=2;
f = figure('color',[1 1 1],'Position',[70.6000 423.4000 1.2584e+03 338.4000]);

subplot(1,2,1)
sub1 =  ezcontourf(Rast,[-1.5 1.5],[-1.5 1.5]);colorbar;set(gca,'fontweight','bold')
hold on

subplot(1,2,2)
sub2 =  ezsurf(Rast,[-1.5 1.5],[-1.5 1.5]);colorbar;set(gca,'fontweight','bold')
colormap jet

%% Input Parameters

[VarMin, VarMax, nVar, VarSize] = SearchSpace();

%% PSOES with default parameters- code snippet (c.3) in the documentation

% 
[sol_c3,sol_c3.ConvHis,paramc3] = psoes(CostFunction,VarSize,VarMin,VarMax);

figure('color',[1 1 1])
ezcontourf(Rast,[-1.5 1.5],[-1.5 1.5])
hold on
colorbar
colormap gray
plot(sol_c3.Position(1),sol_c3.Position(2),'p','markersize',10,'markerfacecolor','g','markeredgecolor','r')
set(gca,'fontweight','bold')


%% PSO with inertia wirght only - code (C.4) in the docs

[sol_c4,sol_c4.ConvHis,paramc4] = psoes(CostFunction,VarSize,VarMin,VarMax,'Mutation',false);


%% PSO with constriction factor

[sol_c5,sol_c5.ConvHis,paramc5] = psoes(CostFunction,VarSize,VarMin,VarMax,'Mutation',false, 'ConstrictionFactor',true);
