function mf=Cost2Dgzz(x,data)
% The cost function performs the forward problem for each potential
% solution, calculates the misfit between the predicted and observed data
% and also applies the regularization constraints, here taken as a simple
% smoothing operator using the gradient function.



data.prism.zdown=x';

gzz = gradGrid3D(data.prism,data.geometry,'gzz','Gamma',-0.5);


mfzz= norm(gzz-data.gzz) ;

S = norm(diff(abs(x)));

mf = mfzz  +0.12* S; 

end