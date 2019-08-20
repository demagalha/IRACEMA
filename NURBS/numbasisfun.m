function B = numbasisfun(i,u,p,U)

n = length(U)-1 -p -1;


B = bsxfun (@(a, b) a+b,i-p, (0:p).').';

end