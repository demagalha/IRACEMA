function connectivity = geo_connectivity(U,pu,qu)

m = length(U)-1;
n = m -pu -1;
ndof = n+1;

nel = size (qu, 2);
nqn = size (qu, 1);


connectivity = zeros (pu+1, nel);
for iel=1:nel
  s = FindSpanList(n, pu, qu(:, iel)', U);
  c = numbasisfun (s, qu(:, iel)', pu, U);
  c = unique(c(:))+1;
  connectivity(1:numel(c), iel) = c;
end

end