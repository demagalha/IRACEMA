function F_local = BuildFLocal2D(R,Jmod,Model,qu,ni,qv,nj,f)
U = Model.U;
pu = Model.pu;
nu = length(U)-pu-1;
u = ((U(ni+1)-U(ni))*qu +U(ni+1)+U(ni))/2;
su = ni;
Nu = BasisFuns(su-1,u,pu,U);

V = Model.V;
pv = Model.pv;
nv = length(V)-pv-1;
v = ((V(nj+1)-V(nj))*qv +V(nj+1)+V(nj))/2;
sv = nj;
Nv = BasisFuns(sv-1,v,pv,V);

Pwx = zeros(pu+1,pv+1);
Pwy = Pwx;
Pwz = Pwx;
Pww = Pwz;
Pw = Model.Pw;

for k=0:pu
    for j=0:pv
        Pwx(k+1,j+1) = Pw(su-pu+k,sv-pv+j).x;
        Pwy(k+1,j+1) = Pw(su-pu+k,sv-pv+j).y;
        Pwz(k+1,j+1) = Pw(su-pu+k,sv-pv+j).z;
        Pww(k+1,j+1) = Pw(su-pu+k,sv-pv+j).w;
    end
end
x = Nu'*Pwx*Nv;
y = Nu'*Pwy*Nv;
z = Nu'*Pwz*Nv;
w = Nu'*Pww*Nv;

F_local = zeros(2,numel(R));
force = f(x);

for aa=1:numel(R)
    F_local(1,aa) = F_local(1,aa) + force(1)*R(aa)*Jmod;
    F_local(2,aa) = F_local(2,aa) + force(2)*R(aa)*Jmod;
end

end