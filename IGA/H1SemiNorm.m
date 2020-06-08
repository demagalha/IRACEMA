function H1_semi = H1SemiNorm(Model)
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
    for i=1:nel
        LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
    end
pu = Model.pu;
pv = Model.pv;
U = Model.U;
V = Model.V;
P = Model.get_point_cell;
 % Get Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu);
    [v, wv] = getGP(pv);
    N_QUAD_U = length(u);
    N_QUAD_V = length(v);
K = zeros(numel(INN(:,1)));
M = K;

N_ELE_DOF = nen;
error = 0;
for e=1:nel
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    if (U(ni+1) == U(ni) || (V(nj+1) == V(nj)))
        continue
    end
    H1_e = 0;
    for i=1:N_QUAD_U
        for j=1:N_QUAD_V
            [R, dR, J] = Shape2D(Model,u(i),v(j),e,P,IEN,INN);
            Jmod = abs(J*wu(i)*wv(j));
            H1_e = H1_e + dot(dR,dR);
        end
    end
    
    % Assemblage
    error = error + H1_e;
end
H1_semi = sqrt(error);
end