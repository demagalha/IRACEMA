function [Sx, Sy, Sz] = SurfacePointRAT(Model,u,v)

n = Model.nu; p = Model.pu; U = Model.U;
m = Model.nv; q = Model.pv; V = Model.V;

%should change this later
Pwx = Model.PX .* Model.weight;
Pwy = Model.PY .* Model.weight;
Pwz = Model.PZ .* Model.weight;
Pww = Model.weight;
%

uspan = FindSpanLinear(n,p,u,U);
Nu = BasisFuns(uspan,u,p,U);
vspan = FindSpanLinear(m,q,v,V);
Nv = BasisFuns(vspan,v,q,V);

Nu = Nu';

PPwx = zeros(p+1,q+1);
PPwy = zeros(p+1,q+1);
PPwz = zeros(p+1,q+1);
PPww = zeros(p+1,q+1);
for k=0:p
    
    for L=0:q
        PPwx(k+1,L+1) = Pwx(uspan-p+k+1,vspan-q+L+1);
        PPwy(k+1,L+1) = Pwy(uspan-p+k+1,vspan-q+L+1);
        PPwz(k+1,L+1) = Pwz(uspan-p+k+1,vspan-q+L+1);
        PPww(k+1,L+1) = Pww(uspan-p+k+1,vspan-q+L+1);
    end
end

w = Nu*PPww*Nv;


Sx = (Nu*PPwx*Nv)/w;
Sy = (Nu*PPwy*Nv)/w;
Sz = (Nu*PPwz*Nv)/w;

end