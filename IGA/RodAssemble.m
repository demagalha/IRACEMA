function [K, M, ID] = RodAssemble(Model,DENSITY,YOUNG_MODULUS)
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
for i = 1:nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),nen,1);
end
%% Assembly
    [x,wx] = getGP(Model.pu); % Quadrature Rule (Gauss Legendre)
    N_QUAD_X = length(x);
    pu = Model.pu;
    rho = DENSITY; % Density
    E = YOUNG_MODULUS; % Elastic 
    N_DOF = numel(INN);
    K = zeros(N_DOF);
    M = K;
    N_ELE_DOF = nen;
    U = Model.U;
    P = Model.get_point_cell;
    

    for e=1:nel % Loop Through Elements
        ni = INN(IEN(1,e),1);
        % Check if element has zero measure
        if U(ni+1) == U(ni)
             continue
        end
        K_e = zeros(nen,nen);
        M_e = K_e;
        for i=1:N_QUAD_X % Loop through quadrature points
            [R, dR, J] = Shape1D(Model,x(i),e,P,IEN,INN);
            Jmod = abs(J*wx(i));
            K_e = K_e + Jmod*E*(dR'*dR);
            M_e = M_e + Jmod*rho*(R'*R);
        end
        idx = LM(:,e)';
        K(idx,idx) = K(idx,idx) + K_e;
        M(idx,idx) = M(idx,idx) + M_e;
    end
end