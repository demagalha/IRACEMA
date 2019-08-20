function iv = FindSpanList(n,p,u,U)

iv = zeros(1,numel(u));

for i=1:numel(u)
    iv(i) = FindSpanLinear(n,p,u(i),U);
end

end
