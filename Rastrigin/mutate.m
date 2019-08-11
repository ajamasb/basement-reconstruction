function xnew = mutate(x,nVar,step,varmin,varmax)
z = zeros(size(x));
for i=1:numel(x)
    z(i)=randn;
end
sigma = step.*(varmax-varmin)./3;
if rand>0.5
    xnew = x + z.*sigma;

else
    xnew = x - z.*sigma;
end
    


end
