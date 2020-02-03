function F_local = BuildFLocal(R,Jmod,Model,qu,ni,qv,nj,qw,nk,f)
U = Model.U;
pu = Model.pu;
u = ((U(ni+1)-U(ni))*qu +U(ni+1)+U(ni))/2;
su = ni;
Nu = BasisFuns(su-1,u,pu,U);

V = Model.V;
pv = Model.pv;
v = ((V(nj+1)-V(nj))*qv +V(nj+1)+V(nj))/2;
sv = nj;
Nv = BasisFuns(sv-1,v,pv,V);

W = Model.W;
pw = Model.pw;
w = ((W(nk+1)-W(nk))*qw +W(nk+1)+W(nk))/2;
sw = nk;
Nw = BasisFuns(sw-1,w,pw,W);
Pw = Model.Pw;
V = CPOINT(0,0,0,0,1);
for L=0:pv
    temp2 = CPOINT(0,0,0,0,1);
    for k=0:pu
        temp = CPOINT(0,0,0,0,1);
        for kk=0:pw
            temp = temp + Nw(kk+1)*Pw(su-pu+k,sv-pv+L,sw-pw+kk);
        end
        temp2 = temp2 + Nu(k+1)*temp;
    end
    V = V + Nv(L+1)*temp2;
end
V = V/V.w;

x = V.x;
y = V.y;
z = V.z;

F_local = zeros(3,numel(R));
force = f(x,y,z);
for aa=1:numel(R)
    F_local(1,aa) = F_local(1,aa) +force(1)*R(aa)*Jmod;
    F_local(2,aa) = F_local(2,aa) +force(2)*R(aa)*Jmod;
    F_local(3,aa) = F_local(3,aa) +force(3)*R(aa)*Jmod;
end
end