function N = AllBasisFuns(span,u,p,U)

N = zeros(p+1,p+1);

for i=0:p
tmp = BasisFuns(span,u,i,U);
for j=1:numel(tmp)
    N(j,i+1) = tmp(j);
end

end