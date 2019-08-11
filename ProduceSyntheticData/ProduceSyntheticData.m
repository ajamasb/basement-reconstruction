% This script produces the synthetic data used in Section 3, Code
% Documentation accompanying papaer Jamasb et al. **

clc;clear;close all
addpath('../Src/')

% first load the model (x,y,z) in meter:
xyz = load('SyntheticModelDepth1.xyz');

dx = mean(diff(xyz(:,1)));

xobs = xyz(:,1)' + dx/2;
yobs = xyz(:,2)';
zobs = zeros(size(xobs)); % perform the calculations at the surface



% Although there is no limitations in defining the model in terms of the
% size of the prisms, for the sake of simplicity, here we design our
% problem to be even-determined (see Menke 1989, ch. 3), which means we
% assign a prism to each observation point. Also, we assume prisms of the
% same lateral extents and densities. Finally, we position the observation
% points exactly at the center of their underlying prisms. 





%% define the structure:

% For simplicity and ease of reproducibility of the resutlts we use a
% 2D synthetic model.

xleft = xyz(:,1);
xright = xleft + dx;


% extend the profiles to avoid artefacts: 
xleft(1) = -10e3;
xright(end) = 10e3;

yleft = -5*max(abs(xyz(:,3))) * ones(size(xleft));
yright = -yleft;

zup = zeros(size(xleft));
zdown = xyz(:,3);

rhoup = -400*ones(size(xleft));
rhodown = rhoup;
gamma = 0.1;


%% calculate the data:
[prism,Obs]=gravGrid3D_makeprism(rhoup,rhodown,xobs,yobs,zobs,xleft,xright,yleft,yright,zup,zdown);

prism.Gamma = gamma;   

% calcualte the gravity effect:
grav = gravGrid3D(prism,Obs,'vector','gamma');


% Calculate Gradiometry tensor:
geometry = prepare_grad_3d(prism,Obs);
grad = gradGrid3D(prism,geometry,'gamma',-0.5);



f1 = figure('color',[1 1 1], 'Position', [521 555 771 294]);
plot(xobs,grav,'k','LineWidth',3);
title('Gravity')
ylabel('mGal')
xlabel('Distance (m)')
set(gca,'FontWeight','bold')
grid minor

f2=figure();
bar(xobs,zdown)




gObs = awgn(grav,100/3,'measured');

figure(1); hold on
plot(xobs,gObs,'LineWidth',3);
legend('True','3% WGN')



gzzObs = awgn(grad.gzz,100/3,'measured');
gxzObs = awgn(grad.gxz,100/3,'measured');
gxxObs = awgn(grad.gxx,100/3,'measured');

f3 = figure('color',[1 1 1 ],'position',  [417 308 1172 223]);

subplot(1,3,1)
plot(xobs,grad.gzz,'k','LineWidth',3);
hold on
plot(xobs,gzzObs,'r','LineWidth',3);
title('G_{zz}')
ylabel('Etvos')
xlabel('Distance (m)')
set(gca,'FontWeight','bold')
grid minor
legend('Synthetic','3% WGN')

subplot(1,3,2)
plot(xobs,grad.gxz,'k','LineWidth',3);
hold on
plot(xobs,gxzObs,'r','LineWidth',3);
title('G_{xz}')
ylabel('Etvos')
xlabel('Distance (m)')
set(gca,'FontWeight','bold')
grid minor
legend('Synthetic','3% WGN')

subplot(1,3,3)
plot(xobs,grad.gxx,'k','LineWidth',3);
hold on
plot(xobs,gxxObs,'r','LineWidth',3);
title('G_{xx}')
ylabel('Etvos')
xlabel('Distance (m)')
set(gca,'FontWeight','bold')
grid minor
legend('Synthetic','3% WGN')


%% save the data

SYNTHETIC.prism = prism;
SYNTHETIC.Obs = Obs;
SYNTHETIC.grav_noNoise = grav;
SYNTHETIC.grad_noNoise = grad;
SYNTHETIC.gObs = gObs;
SYNTHETIC.gzzObs = gzzObs;
SYNTHETIC.gxzObs = gxzObs;
SYNTHETIC.gxxObs = gxxObs;
SYNTHETIC.TrueModel = zdown;



save SYNTHETIC SYNTHETIC




