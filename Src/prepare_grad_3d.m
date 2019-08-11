function geometry = prepare_grad_3d(prism,Obs)
% the function prepares all the constant terms in the gradiometry forward



geometry.sx=size(Obs.xobs);
geometry.sy=size(Obs.yobs);
geometry.sz=size(Obs.zobs);

% if (any(sx~=sy) || any(sx~=sz))
%     warning('The size of observation points are not consistant; Results may not be true!')
% end

RU = prism.RhoUp(:);
RD = prism.RhoDown(:);
XLeft = prism.xleft(:);
XRight = prism.xright(:);
YLeft = prism.yleft(:);
YRight = prism.yright(:);
ZUp = -abs(prism.zup(:)); 
ZDown = -abs(prism.zdown(:));

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


geometry.N = N;
geometry.M = M;

XL = repmat(XLeft,1,N);           % provided that XLeft is a column vector of size Nx1
XR = repmat(XRight,1,N);          % the same
YL = repmat(YLeft,1,N);           % the same
YR = repmat(YRight,1,N);          % the same

geometry.ZU = repmat(ZUp,1,N);             % the same
% ZD = repmat(ZDown,1,N);                    % the same
geometry.RhoUp = repmat(RU,1,N);           % the same
geometry.RhoDown = repmat(RD,1,N);         % the same
% geometry.Gamma = repmat(prism.Gamma(:),M,N);

geometry.XO = repmat(XObs,M,1);            % provided that XObs is a row vector of size 1xM
geometry.YO = repmat(YObs,M,1);            % the same
geometry.ZO = repmat(ZObs,M,1);            % the same


X1 = (XL - geometry.XO)* 1e-3;
X2 = (XR - geometry.XO)* 1e-3;
Y1 = (YL - geometry.YO)* 1e-3;
Y2 = (YR - geometry.YO)* 1e-3;
Z1 = (geometry.ZU - geometry.ZO)* 1e-3;


geometry.X1Y1sq = X1.^2 + Y1.^2;
geometry.X2Y1sq = X2.^2 + Y1.^2;
geometry.X1Y2sq = X1.^2 + Y2.^2;
geometry.X2Y2sq = X2.^2 + Y2.^2;

geometry.X1 = X1; geometry.iX1 = (geometry.X1<1e-10); 
geometry.X2 = X2; geometry.iX2 = (geometry.X2<1e-10);
geometry.Y1 = Y1; geometry.iY1 = (geometry.Y1<1e-10);
geometry.Y2 = Y2; geometry.iY2 = (geometry.Y2<1e-10);
geometry.Z1 = Z1;

geometry.R111=r(X1,Y1,Z1);
geometry.R121=r(X1,Y2,Z1);
geometry.R211=r(X2,Y1,Z1);
geometry.R221=r(X2,Y2,Z1);

a=zeros(size(X1));

geometry.gxz1 = (-gxz1(Y1,geometry.R111) + ...
                        gxz1(Y2,geometry.R121) + ...
                        gxz1(Y1,geometry.R211)  - ...
                        gxz1(Y2,geometry.R221));
                    
geometry.gxz2 = (-gxz2(X1,Y1,Z1,geometry.R111,a)  + ...
                        gxz2(X1,Y2,Z1,geometry.R121,a)  + ...
                        gxz2(X2,Y1,Z1,geometry.R211,a)  - ...
                        gxz2(X2,Y2,Z1,geometry.R221,a) );   
                    
                    
geometry.gzz1 = (-gzz1(X1,Y1,Z1,geometry.R111,a)  + ...
                        gzz1(X1,Y2,Z1,geometry.R121,a)  + ...
                        gzz1(X2,Y1,Z1,geometry.R211,a) - ...
                        gzz1(X2,Y2,Z1,geometry.R221,a) ) ;
                        
                    
                    
                    
geometry.gzz2 =(-gzz2(X1,Y1,geometry.R111,a) + ...
                        gzz2(X1,Y2,geometry.R121,a) + ...
                        gzz2(X2,Y1,geometry.R211,a) - ...
                        gzz2(X2,Y2,geometry.R221,a));                       
                    
                    
                    
geometry.gxx1 = (-gxx1(X1,Y1,Z1,geometry.R111,a) +  ...
                        gxx1(X1,Y2,Z1,geometry.R121,a) + ...
                        gxx1(X2,Y1,Z1,geometry.R211,a) - ...
                        gxx1(X2,Y2,Z1,geometry.R221,a) );...
                        

geometry.gxx2 = (-gxx2(X1,Y1,geometry.R111,a) +...
                        gxx2(X1,Y2,geometry.R121,a) + ...
                        gxx2(X2,Y1,geometry.R211,a) - ...
                        gxx2(X2,Y2,geometry.R221,a));                       %Eotvos   
                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geometry.gyy1 =  (-gyy1(X1,Y1,Z1,geometry.R111,a)  + ...
                        gyy1(X1,Y2,Z1,geometry.R121,a)  + ...
                        gyy1(X2,Y1,Z1,geometry.R211,a)  - ...
                        gyy1(X2,Y2,Z1,geometry.R221,a) );
                        
                    
geometry.gyy2 = (-gyy2(X1,Y1,geometry.R111,a) + ...
                        gyy2(X1,Y2,geometry.R121,a)  + ...
                        gyy2(X2,Y1,geometry.R211,a) - ...
                        gyy2(X2,Y2,geometry.R221,a) );                       %Eotvos  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
geometry.gyz1 =  (-gyz1(X1,geometry.R111) + ...
                        gyz1(X1,geometry.R121) + ...
                        gyz1(X2,geometry.R211) - ...
                        gyz1(X2,geometry.R221) );
                    
                    ...
geometry.gyz2 = (-gyz2(X1,Y1,Z1,geometry.R111,a)+ ...
                        gyz2(X1,Y2,Z1,geometry.R121,a)  + ...
                        gyz2(X2,Y1,Z1,geometry.R211,a)  - ...
                        gyz2(X2,Y2,Z1,geometry.R221,a) );                       %Eotvos 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    

geometry.gxy1 = (-gxy1(Z1,geometry.R111)+ ...
                        gxy1(Z1,geometry.R121)  + ...
                        gxy1(Z1,geometry.R211)  - ...
                        gxy1(Z1,geometry.R221) );
                    ...
geometry.gxy2 =((-geometry.R111)  + (geometry.R121)  + ...
                        (geometry.R211) - (geometry.R221) );                       %Eotvos                      
                                      



end




function R=r(x,y,z)
R = sqrt(x.^2 + y.^2 + z.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gxz1(dy,Rr)
    b =  log( Rr+dy);
    d =b;
end

function d = gxz2(dx,dy,dz,Rr,a)
    c=a;
%     ib = abs(dx)>1e-10;ic = abs(dy)>1e-10; 
%     a(ib)= dx(ib) .* atan( (dz(ib).*dy(ib))./(dx(ib).*Rr(ib)) );
%     c(ic) = dy(ic) .* log( Rr(ic)+dz(ic));
    
    
    a= dx .* atan( (dz.*dy)./(dx.*Rr) );
    c = dy .* log( Rr+dz);
    d =-a+c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gzz1(dx,dy,dz,Rr,a)
    p=a;
%     ia=isnan((dx.*dy)./(dz.*Rr));
%     p(~ia)=  atan( (dx(~ia).*dy(~ia))./(dz(~ia).*Rr(~ia)) );
    
       
    p=  atan( (dx.*dy)./(dz.*Rr) );
    d =p;
end

function d = gzz2(dx,dy,Rr,a)
    b=a;c=a;
%     ib = abs(dx)>1e-10;ic = abs(dy)>1e-10;    
%     b(ic) = dy(ic) .* log( Rr(ic)+dx(ic));
%     c(ib) = dx(ib) .* log( Rr(ib)+dy(ib));
    
   
    b = dy .* log( Rr+dx);
    c = dx.* log( Rr+dy);
    d =b+c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gxx1(dx,dy,dz,Rr,a)
%     p=a;
%     ia=isnan((dz.*dy)./(dx.*Rr));
%     p(~ia)=  atan( (dz(~ia).*dy(~ia))./(dx(~ia).*Rr(~ia)) );
%     d =p;

   p=a;
   
    p=  atan( (dz.*dy)./(dx.*Rr) );
    d =p;
end

function d = gxx2(dx,dy,Rr,a)
%     c=a;
%     ib = abs(dx)>1e-10;
%     c(ib) = dx(ib) .* log( Rr(ib)+dy(ib));
%     d =c;

    c=a;
    
    c = dx .* log( Rr+dy);
    d =c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gyy1(dx,dy,dz,Rr,a)
%     p=a;
%     ia=isnan((dz.*dx)./(dy.*Rr));
%     p(~ia)=  atan( (dz(~ia).*dx(~ia))./(dy(~ia).*Rr(~ia)) );
%     d =p;

    p=a;
    
    p=  atan( (dz.*dx)./(dy.*Rr) );
    d =p;
end

function d = gyy2(dx,dy,Rr,a)
%     c=a;
%     ib = abs(dy)>1e-10;
%     c(ib) = dy(ib) .* log( Rr(ib)+dx(ib));
%     d =c;

    c= dy .* log( Rr+dx);
    d =c;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gyz1(dx,Rr)
    b =  log(Rr+dx);
    d =b; 
end

function d = gyz2(dx,dy,dz,Rr,a)
%     c=a;
%     ib = abs(dx)>1e-10;ic = abs(dy)>1e-10; 
%     a(ic)= dy(ic) .* atan( (dx(ic).*dz(ic))./(dy(ic).*Rr(ic)) );
%     c(ib) = dx(ib) .* log( Rr(ib)+dz(ib));
%     d =-a+c;
    p=a;

    p=  atan( (dz.*dx)./(dy.*Rr) );
    d =p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gxy1(dz,Rr)
    b =   log(Rr+dz);
    d =b; 
end
