function [ VarMin, VarMax, nVar, VarSize] = SearchSpace()
% search space for the minimization of the Rastrigin Function.
% see the docs.
%
% 
% Ali Jamasb,
% Institute of Geophysics, University of Tehran, Iran.
% ajamasb@ut.ac.ir
% Jul. 14, 2019
%

nVar = 2;  % the number of model parameters i.e. (x & y)

fprintf(' The Number of Model Parameters: %d\n',nVar)

VarSize = [1,nVar];


VarMin = repmat(-10,VarSize);  % the lower limits of the search space
VarMax = repmat(10,VarSize);    % the upper limits of the search space

end

