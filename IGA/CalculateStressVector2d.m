function stress = CalculateStressVector2d(P,Model,DX,DY,u,v,YOUNG,POISSON )

nu = Model.nu;
nv = Model.nv;
pu = Model.pu;
pv = Model.pv;
U = Model.U;
V = Model.V;

dR = CalculateDers2D(Model,u,v,P);
% if any(isnan([dR(:,1);dR(:,2)])) || any(isinf([dR(:,1);dR(:,2)]))
%     if u ~= 1 && v ~= 1
%         dR = CalculateDers2D(Model,u+0.01,

uspan = FindSpanLinear(nu,pu,u,U); vspan = FindSpanLinear(nv,pv,v,V);



del_u = zeros(2,2);
loc = 0;
for j=0:pv
    for i=0:pu
        loc = loc + 1;
        del_u(1,1) = del_u(1,1) + dR(loc,1)*DX(uspan-pu+i+1,vspan-pv+j+1);
        del_u(1,2) = del_u(1,2) + dR(loc,2)*DX(uspan-pu+i+1,vspan-pv+j+1);
        
        del_u(2,1) = del_u(2,1) + dR(loc,1)*DY(uspan-pu+i+1,vspan-pv+j+1);
        del_u(2,2) = del_u(2,2) + dR(loc,2)*DY(uspan-pu+i+1,vspan-pv+j+1);
    end
end

strain_vector = [del_u(1,1); del_u(2,2); del_u(1,2) + del_u(2,1)];
C = (YOUNG/(1-POISSON^2))*[1 POISSON 0; POISSON 1 0; 0 0 0.5*(1-POISSON)];

stress = C*strain_vector;

end
        