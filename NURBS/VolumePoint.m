function V = VolumePoint(nu,pu,U,nv,pv,V,nw,pw,W,Pw,u,v,w)

%input: nu,pu,U <--- the u direction, nu being n+1 control points in the
%direction u, pu the degree in that direction and U the knot vector
%and so on for the other two directions... Pw is the "3D matrix" of cpts in the homogenous space
%u,v,w parameters for the point to be valued at
%
%output: V, the point on euclidian space

uspan = FindSpanLinear(nu,pu,u,U);
Nu = BasisFuns(uspan,u,pu,U);
vspan = FindSpanLinear(nv,pv,v,V);
Nv = BasisFuns(vspan,v,pv,V);
wspan = FindSpanLinear(nw,pw,w,W);
Nw = BasisFuns(wspan,w,pw,W);

V = CPOINT(0,0,0,0,1); 
for L=0:pv
    temp2 = CPOINT(0,0,0,0,1);
    for k=0:pu
        temp = CPOINT(0,0,0,0,1);
        for kk=0:pw
            temp = temp + Nw(kk+1)*Pw(uspan-pu+k+1,vspan-pv+L+1,wspan-pw+kk+1);
        end
        temp2 = temp2 + Nu(k+1)*temp;
    end
    V = V + Nv(L+1)*temp2;
end
V = V/V.w;
end