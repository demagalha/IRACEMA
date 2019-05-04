function [k,s] = FindSpanMult(n,p,u,U)

k = FindSpanLinear(n,p,u,U);
s = Mult(n,p,u,U);

end