function g=gravGrid3D(prism,Obs,OutMode,FixedMode)
%
% function g=grav_vldg(RhoUp,RhoDown,XObs,YObs,ZObs,Xp,Yp,Zp) calculates the gravity effect of
% a right rectangular prism whose density varies linearly in vertical
% direction. The formula is from  Gallardo-Delgado et al., Geophysics, 68
% (2002), p.949.
%
% ========================================================================
% Input  :
% * * * * * *
% prism: a structure array containing the densities and the coordinates of
%         the prisms;
% Obs:   a structure array containing the coordinates of the observation
%        points.
% * * * * *  * * * * ****************************** * * * * * * * * * * * *
% Produce the inputs using gravGrid3D_makeprism code - if not available copy
% and paste the following in a seperate m file and save:
% % --------------------------------------------------------------------------
% function [prism,Obs]=gravGrid3D_makeprism(ru,rd,xo,yo,zo,xl,xr,yl,yr,zu,zd)
% 
% prism.RhoUp = ru;
% prism.RhoDown =rd ;
% prism.xleft = xl;
% prism.xright = xr;
% prism.yleft = yl;
% prism.yright = yr;
% prism.zup = zu;
% prism.zdown =  zd; 
% 
% Obs.xobs = xo ;
% Obs.yobs = yo;
% Obs.zobs = zo;
% end
% --------------------------------------------------------------------------

%
% OutMode: a string: could be 'row', 'vector', or 'original'
%          1)'row'/'vector': in this case the output will be a Nx3 array:
%          [x,y,gi]
%           2)'original': useful for grids: gives the output in the same
%           format that the observation points are given; if the x and y
%           coordinates of the observation points are produced with
%           meshgrid, then g will be of the same size;


% ***********************************************************************
% **(could be a point, line or a grid -- for a line (profile) the other**
% **horizotal direction must be extended to infinity)********************
% ***********************************************************************
%

%
% * Gamma: the slope of the increasing density (RHO(z)=Gamma*z + RHO1)
% * or ->: Rho_up and Rho_down in (Kg/m3)
% ========================================================================
% --->> ALL the input distances must be in meter
%
% Output:
% * g   : the gravity in mGal


validatestring(OutMode,{'row','vector','original'});
validatestring(FixedMode,{'density','gamma'});

sx=size(Obs.xobs);
sy=size(Obs.yobs);
sz=size(Obs.zobs);

if (any(sx~=sy) || any(sx~=sz))
    warning('The size of observation points are not consistant; Results may not be true!')
end

RU = prism.RhoUp(:);
RD = prism.RhoDown(:);
XLeft = prism.xleft(:);
XRight = prism.xright(:);
YLeft = prism.yleft(:);
YRight = prism.yright(:);
ZUp = prism.zup(:);
ZDown = prism.zdown(:);

XObs = Obs.xobs(:)';
YObs = Obs.yobs(:)';
ZObs = Obs.zobs(:)';
% *********** * * * * * * * * * * * * * * * * * * * * * * * * *************
% ASSUMING THAT EVERYTHING IS SORTED RIGHT (SIZE WISE) THEN WE SHOULD PROCEED
% *********** * * * * * * * * * * * * * * * * * * * * * * * * *************


% IF WE HAVE N OBSERVATION POINTS AND M PRISMS:
% THEN WE NEED EVERYTHING TO BE OF SIZE MxN
N = numel(XObs);
M = numel(XLeft);

XL = repmat(XLeft,1,N);           % provided that XLeft is a column vector of size Nx1
XR = repmat(XRight,1,N);          % the same
YL = repmat(YLeft,1,N);           % the same
YR = repmat(YRight,1,N);          % the same
ZU = repmat(ZUp,1,N);             % the same
ZD = repmat(ZDown,1,N);           % the same
RhoUp = repmat(RU,1,N);           % the same
RhoDown = repmat(RD,1,N);         % the same

XO = repmat(XObs,M,1);            % provided that XObs is a row vector of size 1xM
YO = repmat(YObs,M,1);            % the same
ZO = repmat(ZObs,M,1);            % the same


X1 = (XL - XO) * 1e-3;
X2 = (XR - XO)* 1e-3;
Y1 = (YL - YO)* 1e-3;
Y2 = (YR - YO)* 1e-3;
Z1 = (ZU - ZO)* 1e-3;
Z2 = (ZD - ZO)* 1e-3;















% 
% if (Zp(2) == Zp(1))
%     g = 0;
%     warning(' The Zp coordinates are the same, gravity is returned zero')
%     stop
% end

switch FixedMode
    case 'density'
        Gamma = (RhoDown - RhoUp)./(ZD - ZU);
    case 'gamma'
        Gamma = repmat(prism.Gamma(:),M,N);
        
end

Rho0  = (RhoUp-Gamma.*ZU)*1e-3;
RhoZ0 = (Rho0+Gamma.*ZO*1e-3 );



a=zeros(size(X1));

format long
%% Constants:
G = 6.6732e-11;

g = (RhoZ0.*G.*1e11) .* (-gxyz(X1,Y1,Z1,r(X1,Y1,Z1),a) + gxyz(X1,Y1,Z2,r(X1,Y1,Z2),a) + ...
                        gxyz(X1,Y2,Z1,r(X1,Y2,Z1),a) - gxyz(X1,Y2,Z2,r(X1,Y2,Z2),a) + ...
                        gxyz(X2,Y1,Z1,r(X2,Y1,Z1),a) - gxyz(X2,Y1,Z2,r(X2,Y1,Z2),a) - ...
                        gxyz(X2,Y2,Z1,r(X2,Y2,Z1),a) + gxyz(X2,Y2,Z2,r(X2,Y2,Z2),a)) ...
                        + (Gamma.*G.*1e11).*(-gxyzGamma(X1,Y1,Z1,r(X1,Y1,Z1),a) +...
                        gxyzGamma(X1,Y1,Z2,r(X1,Y1,Z2),a) + ...
                        gxyzGamma(X1,Y2,Z1,r(X1,Y2,Z1),a) -...
                        gxyzGamma(X1,Y2,Z2,r(X1,Y2,Z2),a) + ...
                        gxyzGamma(X2,Y1,Z1,r(X2,Y1,Z1),a) -...
                        gxyzGamma(X2,Y1,Z2,r(X2,Y1,Z2),a) - ...
                        gxyzGamma(X2,Y2,Z1,r(X2,Y2,Z1),a) +...
                        gxyzGamma(X2,Y2,Z2,r(X2,Y2,Z2),a));

                    
                    
                    
                    
                    
                    switch OutMode
                        case {'row','vector'}
                            g = [sum(g)];
                        case 'original'
                            g = reshape(sum(g),sx);
                    end
format short
end

function R=r(x,y,z)
R = sqrt(x.^2 + y.^2 + z.^2);
end

function d = gxyz(dx1,dy1,dz1,R111,a)
b=a;c=a;

a= dz1 .* atan( (dx1.*dy1)./(dz1.*R111) );
b = dx1 .* log( R111+dy1 );
c= dy1 .* log( R111+dx1);
d =a- b- c;
end
function d = gxyzGamma(dx1,dy1,dz1,R111,a)
% b=a;c=a;l=a;
% ia = abs(dz1)>1e-10; ib= abs(dx1)>1e-10 ; ic =abs(dy1)>1e-10;il =(abs(dx1)>1e-10 & abs(dy1)>1e-10);
a=0.5*dz1.^2 .* atan( (dx1.*dy1)./(dz1.*R111) );
b=0.5*dx1.^2 .* atan( (dz1.*dy1)./(dx1.*R111) );
c=0.5*dy1.^2 .* atan( (dx1.*dz1)./(dy1.*R111) );
l = dx1.*dy1 .* log( R111+dz1);
d =a-b-c+ l;

end
