function S = SurfacePointRAT3(n,p,U,m,q,V,Pw,u,v,S)

uspan = FindSpanLinear(n,p,u,U);
Nu = BasisFuns(uspan,u,p,U);

vspan = FindSpanLinear(m,q,v,V);
Nv = BasisFuns(vspan,v,q,V);

Nu = Nu';

PPw.x = zeros(p+1,q+1);
PPw.y = zeros(p+1,q+1);
PPw.z = zeros(p+1,q+1);
PPw.w = zeros(p+1,q+1);
for k=0:p
    
    for L=0:q
        PPw.x(k+1,L+1) = Pw(uspan-p+k+1,vspan-q+L+1).x;
        PPw.y(k+1,L+1) = Pw(uspan-p+k+1,vspan-q+L+1).y;
        PPw.z(k+1,L+1) = Pw(uspan-p+k+1,vspan-q+L+1).z;
        PPw.w(k+1,L+1) = Pw(uspan-p+k+1,vspan-q+L+1).w;
    end
end

w = Nu*PPw.w*Nv;


S.x = (Nu*PPw.x*Nv)/w;
S.y = (Nu*PPw.y*Nv)/w;
S.z = (Nu*PPw.z*Nv)/w;

end