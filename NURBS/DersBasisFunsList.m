function dB = DersBasisFunsList(i,u,p,n,U)

dB = cell(numel(i),1);

for ii=1:numel(i)
    dB{ii} = DersBasisFun(i(ii),u(ii),p,n,U);
end

end
