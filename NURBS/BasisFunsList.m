function B = BasisFunsList(i,u,p,U)

B = zeros(numel(u),p+1);

for ii=1:numel(u)
    B(ii,:) = BasisFuns(i(ii),u(ii),p,U);
end

end
