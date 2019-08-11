function [prism,Obs]=gravGrid3D_makeprism(ru,rd,xo,yo,zo,xl,xr,yl,yr,zu,zd)

prism.RhoUp = ru;
prism.RhoDown =rd ;
prism.xleft = xl;
prism.xright = xr;
prism.yleft = yl;
prism.yright = yr;
prism.zup = zu;
prism.zdown =  zd; 

Obs.xobs = xo ;
Obs.yobs = yo;
Obs.zobs = zo;
end
