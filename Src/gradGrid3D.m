function [out]=gradGrid3D(prism,geometry,varargin)
% function out=gradGrid3D(prism,geometry,varargin) calculates the gravity
% gradiometry tensor for either a DEM or a basin. The function uses right
% rectangular prisms with a linear density gradient to discretize the mass 
% for which the data is to be calculated.
% See the accompanying paper (Jamasb et al. ***, ***) for formulas.
%
%
% INPUT PARAMETERS (All Units are in SI)-> Distances in (m); density(kg/m3)
%-------------------------------------------------------------------------
% * prism(Required)      :  a struct holding the coordinates of the prisms
% * geometry(Required)   :  a struct holding the lateral
% 
% * Output (Optional)    :  Default output is a struct with six feilds for each
%                           gravity gradiometry component - otherwise choose
%                           between:
%                           {'gxx','gyy','gzz','gxz','gzy','gxy','all'} 
% 
% * 'Gamma',value        : The linear gradient can be defined inderectly
%  (Optional Pair)        from the top and bottom densities (default) or
%                          it can be given directly. 
% 
% ------------------------------------------------------------------------
% Please note that the order of input follows MATLAB principles:
% First enter the Required inputs, then optionals, and finally optional pairs.
% ------------------------------------------------------------------------
% 
% Example: 
% out=gradGrid3D(prism,geometry);            calculates the full tensor 
% 
% out=gradGrid3D(prism,geometry,'gxx');      calculates only gxx component
% 
% out=gradGrid3D(prism,geometry,'Gamma',0);  calculates the full tensor for
%                                            a prism with constant density (i.e. Gamma=0)
%
% out=gradGrid3D(prism,geometry,'gyz','Gamma',0.02); calculates the gyz component of for a prism with a
%                                                    linear density gradient of 0.02
%
%
%--------------------------------------------------------------------------
% Produce the inputs using the function:
% "[prism,obs]=gravGrid3D_makeprism(RhoUp,RhoDown,xobs,yobs,zobs,xleft,xright,yleft,yright,zup,zdown)" 
%  and "prepare_grad_3d", of which only the first needs user-defined
%  inputs including (in case of M prisms and N observation points):
% 
%  1) The densities at the top and bottom of prisms: RhoUp & RhoDown Mx1 vectors
% 
%  2) The coordinates of the observation points: xobs, yobs,zobs  1xN vector
%    
%  3) The lateral extents of the prisms: xleft, xright, yleft, yright   Mx1
%  vectors
% 
%  4) and the upper and lower limits of the prisms: zup, zdown Mx1 vectors
% 
% 
% 
% 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Please note that the code gradGrid3D is primarily designed for interface
% (i.e. Basement, Moho, etc.) reconstruction. Such a problem is non-linear
% because the unknown of the problem is z which has a non-linear
% relashionship with gravity. Therefore, such inversion shall be solved
% using iterative methods either deterministic or stochastic. In the
% former, the forward nedds to be calculated at least twice in each
% iteration (including the calculation of Freshet Matrix) while in
% stochastics, the forward would be calculated many times in each
% iteration. Thus, it is beneficial to perform that part of the forward
% which is fixed during the inversion (those terms that do not include the
% depth at the bottom) beforehand. This would reduce the calculation time
% by half. That part here is calculated using the function geometry = prepare_grad_3d(prism,Obs)
% the result of which is given as input to this code.
% For calculating the forwards only once, for instance in case of a DEM,
% use out=gradGrid3DDEN(prism,Obs).
% 

% Developed by:
% A. Jamasb (ajamasb@ut.ac.ir)
% Institute of Geophysics,
% University of Tehran,
% Iran.
 





narginchk(2,5)



p = inputParser;

DefaultOutput = 'all';
PossibleOutputs = {'gxx','gyy','gzz','gxz','gzx','gyz','gzy','gxy','gyx','all'};

addRequired(p,'prism', @(x)(validateattributes(x,{'struct'},{'nonempty'})));
addRequired(p,'geometry',@(x)(validateattributes(x,{'struct'},{'nonempty'})));
addOptional(p, 'Output',DefaultOutput,@(x)(any(validatestring(x,PossibleOutputs))));
addParameter(p,'Gamma',[],@(x)(validateattributes(x,{'double'},{'nonempty'})));
parse(p,prism,geometry,varargin{:})




ZD = -abs(prism.zdown(:));
Z2 = (repmat(ZD,1,geometry.N) - geometry.ZO)* 1e-3; 
Z2sq = Z2.^2;
        

if isempty(p.Results.Gamma)
    
    
%     geometry.Gamma = ((geometry.RhoDown - geometry.RhoUp)./((Z2 - geometry.Z1)))*1e-3;
      geometry.Gamma = (( geometry.RhoDown -  geometry.RhoUp)./(Z2 -  geometry.Z1))*1e-3;
    
else
   
    geometry.Gamma = repmat(p.Results.Gamma,geometry.M,geometry.N);
     
end
   
        
Rho0  = (geometry.RhoUp-(geometry.Gamma).*geometry.ZU)*1e-3;
RhoZ0 = (Rho0+(geometry.Gamma).*geometry.ZO*1e-3 );

a=zeros(size(geometry.X1));

R112 = sqrt(geometry.X1Y1sq + Z2sq);

R122 = sqrt(geometry.X1Y2sq + Z2sq);

R212= sqrt(geometry.X2Y1sq + Z2sq);

R222= sqrt(geometry.X2Y2sq + Z2sq);





G = 6.6732e-11;




switch p.Results.Output
    
    case {'gxz','gzx',}

gxz = -(RhoZ0.*G.*1e12) .* ( + gxz1(geometry.Y1,R112)  ...
                         - gxz1(geometry.Y2,R122)  ...
                         - gxz1(geometry.Y1,R212)  ...
                         + gxz1(geometry.Y2,R222) + geometry.gxz1) ...
                        + (geometry.Gamma.*G.*1e12).*( + gxz2(geometry.X1,geometry.Y1,Z2,R112,a) ...
                         - gxz2(geometry.X1,geometry.Y2,Z2,R122,a)  ...
                         - gxz2(geometry.X2,geometry.Y1,Z2,R212,a) ...
                         + gxz2(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gxz2);                       %Eotvos
                                      

finalxx = sum(gxz);
out = reshape(finalxx,geometry.sx);


 case {'gzz'}
gzz = (RhoZ0.*G.*1e12) .* ( gzz1(geometry.X1,geometry.Y1,Z2,R112,a) ...
                         - gzz1(geometry.X1,geometry.Y2,Z2,R122,a)  ...
                         - gzz1(geometry.X2,geometry.Y1,Z2,R212,a) ...
                        + gzz1(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gzz1) ...
                        +  (geometry.Gamma.*G.*1e12).*(gzz2(geometry.X1,geometry.Y1,R112,a)...
                         - gzz2(geometry.X1,geometry.Y2,R122,a) ...
                        - gzz2(geometry.X2,geometry.Y1,R212,a) ...
                         + gzz2(geometry.X2,geometry.Y2,R222,a) + geometry.gzz2);                       %Eotvos
                    
finalzz = sum(gzz);

out = reshape((finalzz),geometry.sx);


 case {'gxx'}
                    
gxx = (RhoZ0.*G.*1e12) .* (  gxx1(geometry.X1,geometry.Y1,Z2,R112,a) + ...
                        - gxx1(geometry.X1,geometry.Y2,Z2,R122,a) + ...
                         - gxx1(geometry.X2,geometry.Y1,Z2,R212,a) - ...
                         + gxx1(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gxx1) ...
                        - (geometry.Gamma.*G.*1e12).*(  gxx2(geometry.X1,geometry.Y1,R112,a)  ...
                         - gxx2(geometry.X1,geometry.Y2,R122,a)  ...
                         - gxx2(geometry.X2,geometry.Y1,R212,a) ...
                        + gxx2(geometry.X2,geometry.Y2,R222,a)+ geometry.gxx2);   		 						%Eotvos 

                    finalxx = sum(gxx);

out = reshape((finalxx),geometry.sx);
                    
 case {'gyy'}                   
                    
gyy = (RhoZ0.*G.*1e12) .* (  gyy1(geometry.X1,geometry.Y1,Z2,R112,a)  ...
                        - gyy1(geometry.X1,geometry.Y2,Z2,R122,a) + ...
                         - gyy1(geometry.X2,geometry.Y1,Z2,R212,a)  ...
                         + gyy1(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gyy1) ...
                        - (geometry.Gamma.*G.*1e12).*( gyy2(geometry.X1,geometry.Y1,R112,a)  ...
                       - gyy2(geometry.X1,geometry.Y2,R122,a)  ...
                         - gyy2(geometry.X2,geometry.Y1,R212,a)  ...
                         + gyy2(geometry.X2,geometry.Y2,R222,a) + geometry.gyy2);  
                     
                    finalyy = sum(gyy);

out = reshape((finalyy),geometry.sx);                     
    
                     
  case {'gyz','gzy'}                     
                     
gyz = -(RhoZ0.*G.*1e12) .* ( gyz1(geometry.X1,R112)  ...
                        - gyz1(geometry.X1,R122)  ...
                        - gyz1(geometry.X2,R212)  ...
                         + gyz1(geometry.X2,R222) + geometry.gyz1) ...
                        + (geometry.Gamma.*G.*1e12).*(gyz2(geometry.X1,geometry.Y1,Z2,R112,a)  ...
                        - gyz2(geometry.X1,geometry.Y2,Z2,R122,a)  ...
                        - gyz2(geometry.X2,geometry.Y1,Z2,R212,a)  ...
                        + gyz2(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gyz2);                       %Eotvos        
   

                    finalyz = sum(gyz);

out = reshape((finalyz),geometry.sx);                      
                    
 case {'gyx','gxy'}                     
gxy = -(RhoZ0.*G.*1e12) .* ( gxy1(Z2,R112) ...
                         - gxy1(Z2,R122)  ...
                        - gxy1(Z2,R212)  ...
                         + gxy1(Z2,R222) + + geometry.gxy1) ...
                        - (geometry.Gamma.*G.*1e12).*(  (R112) - (R122)  - (R212)  + (R222) + geometry.gxy2); 
                    
                    finalxy = sum(gxy);

out = reshape((finalxy),geometry.sx);                      
                    
                    
                    
    case {'all'}
  
   

gxz = -(RhoZ0.*G.*1e12) .* ( + gxz1(geometry.Y1,R112)  ...
                         - gxz1(geometry.Y2,R122)  ...
                         - gxz1(geometry.Y1,R212)  ...
                         + gxz1(geometry.Y2,R222) + geometry.gxz1) ...
                        + (geometry.Gamma.*G.*1e12).*( + gxz2(geometry.X1,geometry.Y1,Z2,R112,a) ...
                         - gxz2(geometry.X1,geometry.Y2,Z2,R122,a)  ...
                         - gxz2(geometry.X2,geometry.Y1,Z2,R212,a) ...
                         + gxz2(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gxz2);                       %Eotvos
                                      

finalxx = sum(gxz);
gxz = reshape(finalxx,geometry.sx);



gzz = (RhoZ0.*G.*1e12) .* ( gzz1(geometry.X1,geometry.Y1,Z2,R112,a) ...
                         - gzz1(geometry.X1,geometry.Y2,Z2,R122,a)  ...
                         - gzz1(geometry.X2,geometry.Y1,Z2,R212,a) ...
                        + gzz1(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gzz1) ...
                        +  (geometry.Gamma.*G.*1e12).*(gzz2(geometry.X1,geometry.Y1,R112,a)...
                         - gzz2(geometry.X1,geometry.Y2,R122,a) ...
                        - gzz2(geometry.X2,geometry.Y1,R212,a) ...
                         + gzz2(geometry.X2,geometry.Y2,R222,a) + geometry.gzz2);                       %Eotvos
                    
finalzz = sum(gzz);

gzz = reshape((finalzz),geometry.sx);



                    
gxx = (RhoZ0.*G.*1e12) .* (  gxx1(geometry.X1,geometry.Y1,Z2,R112,a) + ...
                        - gxx1(geometry.X1,geometry.Y2,Z2,R122,a) + ...
                         - gxx1(geometry.X2,geometry.Y1,Z2,R212,a) - ...
                         + gxx1(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gxx1) ...
                        - (geometry.Gamma.*G.*1e12).*(  gxx2(geometry.X1,geometry.Y1,R112,a)  ...
                         - gxx2(geometry.X1,geometry.Y2,R122,a)  ...
                         - gxx2(geometry.X2,geometry.Y1,R212,a) ...
                        + gxx2(geometry.X2,geometry.Y2,R222,a)+ geometry.gxx2);   		

 finalxx = sum(gxx);

gxx = reshape((finalxx),geometry.sx);                   
                 
                    
gyy = (RhoZ0.*G.*1e12) .* (  gyy1(geometry.X1,geometry.Y1,Z2,R112,a)  ...
                        - gyy1(geometry.X1,geometry.Y2,Z2,R122,a) + ...
                         - gyy1(geometry.X2,geometry.Y1,Z2,R212,a)  ...
                         + gyy1(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gyy1) ...
                        - (geometry.Gamma.*G.*1e12).*( gyy2(geometry.X1,geometry.Y1,R112,a)  ...
                       - gyy2(geometry.X1,geometry.Y2,R122,a)  ...
                         - gyy2(geometry.X2,geometry.Y1,R212,a)  ...
                         + gyy2(geometry.X2,geometry.Y2,R222,a) + geometry.gyy2);  
    
                     
 finalyy = sum(gyy);

gyy = reshape((finalyy),geometry.sx);                   
                     
gyz = -(RhoZ0.*G.*1e12) .* ( gyz1(geometry.X1,R112)  ...
                        - gyz1(geometry.X1,R122)  ...
                        - gyz1(geometry.X2,R212)  ...
                         + gyz1(geometry.X2,R222) + geometry.gyz1) ...
                        + (geometry.Gamma.*G.*1e12).*(gyz2(geometry.X1,geometry.Y1,Z2,R112,a)  ...
                        - gyz2(geometry.X1,geometry.Y2,Z2,R122,a)  ...
                        - gyz2(geometry.X2,geometry.Y1,Z2,R212,a)  ...
                        + gyz2(geometry.X2,geometry.Y2,Z2,R222,a) + geometry.gyz2);                       %Eotvos        
   
        finalyz = sum(gyz);

gyz= reshape((finalyz),geometry.sx);            
                  
gxy = -(RhoZ0.*G.*1e12) .* ( gxy1(Z2,R112) ...
                         - gxy1(Z2,R122)  ...
                        - gxy1(Z2,R212)  ...
                         + gxy1(Z2,R222) + + geometry.gxy1) ...
                        - (geometry.Gamma.*G.*1e12).*(  (R112) - (R122)  - (R212)  + (R222) + geometry.gxy2);                     
                     
                     
 
 finalxy = sum(gxy);

gxy = reshape((finalxy),geometry.sx);    


    out.gxx = gxx;
    out.gxz = gxz;
    out.gxy = gxy;
    out.gzz = gzz;
    out.gyy = gyy;
    out.gyz = gyz;
    
                     
end   



    
    
    
    
    


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%