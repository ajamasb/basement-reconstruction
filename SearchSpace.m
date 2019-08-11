function [ VarMin, VarMax, nVar, VarSize, data] = SearchSpace(input)

% Some of the operations here are trivial in the case of this synthetic
% example. However, the following procedure is suitable for real data.

data.prism = input.prism;
data.obs = input.Obs;
data.gzz = input.gzzObs;
data.g = input.gObs;

nVar = numel(input.prism.xleft);
VarSize = [1,nVar];
VarMin = zeros([1, nVar]);
VarMax = repmat(700,1,nVar);

geometry = prepare_grad_3d(input.prism,input.Obs);
data.geometry =geometry; 

end

