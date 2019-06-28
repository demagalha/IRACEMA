function SKL = SurfDerivAlg1(n,p,U,m,q,V,P,u,v,d)

du = min(d,p);
for k=p+1:d
    for L=0:d-k
        SKL(k+1,L+1) = 0;
    end
end

dv = min(d,q);
for L=q+1:d
    for k=0:d-L
        SKL(k+1,L+1) = 0;
    end
end

uspan = FindSpanLinear(n,p,u,U);
Nu = DersBasisFun(uspan,u,p,du,U);

vspan = FindSpanLinear(m,q,v,V);
Nv = DersBasisFun(vspan,v,q,dv,V);

for k=0:du
    for s=0:q
        temp(s+1) = 0;
        for r=0:p
            temp(s+1) = temp(s+1) + Nu(k+1,r+1)*P(uspan-p+r+1,vspan-q+s+1);
        end
    end
    dd = min(d-k,dv);
    for L=0:dd
        SKL(k+1,L+1) = 0;
        for s=0:q
            SKL(k+1,L+1) = SKL(k+1,L+1) + Nv(L+1,s+1)*temp(s+1);
        end
    end
end
end