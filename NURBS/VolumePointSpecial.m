function [Vx,Vy,Vz] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u,v,w)

nu = Model.nu; pu = Model.pu; U = Model.U;
nv = Model.nv; pv = Model.pv; V = Model.V;
nw = Model.nw; pw = Model.pw; W = Model.W; 

uspan = FindSpanLinear(nu,pu,u,U);
Nu = BasisFuns(uspan,u,pu,U);
vspan = FindSpanLinear(nv,pv,v,V);
Nv = BasisFuns(vspan,v,pv,V);
wspan = FindSpanLinear(nw,pw,w,W);
Nw = BasisFuns(wspan,w,pw,W);

Pww = Model.weight;

V = [0,0,0,0];
for L=0:pv
    temp2 = [0,0,0,0];
    for k=0:pu
        temp = [0,0,0,0];
        for kk=0:pw
            Pw = [Pwx(uspan-pu+k+1,vspan-pv+L+1,wspan-pw+kk+1), Pwy(uspan-pu+k+1,vspan-pv+L+1,wspan-pw+kk+1), Pwz(uspan-pu+k+1,vspan-pv+L+1,wspan-pw+kk+1), Pww(uspan-pu+k+1,vspan-pv+L+1,wspan-pw+kk+1)];
            temp = temp + Nw(kk+1)*Pw;
        end
        temp2 = temp2 + Nu(k+1)*temp;
    end
    V = V + Nv(L+1)*temp2;
end
V = V/V(4);
Vx = V(1); Vy = V(2); Vz = V(3);
end
