function F = NeumannBoundaryCondition(Model,F,BoundaryElements)
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
    for i=1:nel
        LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
    end
pu = Model.pu;
pv = Model.pv;
p = [pu; pv];
U = Model.U;
V = Model.V;
Knots = {U; V};
P = Model.get_point_cell;
 % Get Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu);
    [v, wv] = getGP(pv);
    quadpoints = {u, wu; v, wv};
    N_QUAD_U = length(u);
    N_QUAD_V = length(v);
    N_QUAD = [N_QUAD_U; N_QUAD_V];
N_ELE_DOF = nen;
for ee=1:length(BoundaryElements)
    e = BoundaryElements(ee,1);
    ni = INN(IEN(1,e),:);
    parametric_direction = BoundaryElements(ee,2); % Direction of the Boundary
    direction = setdiff([1 2], parametric_direction); % Direction of Integration
    pidx = ni(parametric_direction)-p(parametric_direction):ni(parametric_direction);
    % This is done so we integrate the boundary in the OTHER direction
        q = quadpoints{direction,1};
        w = quadpoints{direction,2};
        Domain = Knots{direction};
        pD = p(direction);
        order = p(direction);
        n = ni(direction);
        h = Domain(n+1) - Domain(n); % Element size-scale
        % Check if element has zero measure
        if h < sqrt(eps) % if element has zero measure
            continue
        end
    bound_val = BoundaryElements(ee,3);
    lift = BoundaryElements(ee,4);
    F_e = zeros(N_ELE_DOF,1);
    K_e = zeros(N_ELE_DOF);
    for i=1:N_QUAD(direction)        
        if direction == 1
            qu = q(i);
            qv =(BoundaryElements(ee,3)*2 -1);
        elseif direction == 2
            qu = (BoundaryElements(ee,3)*2 -1);
            qv = q(i);
        end
        [R, ~, J] = Shape2D(Model,qu,qv,e,P,IEN,INN);
        F_e = F_e + abs(J*w(i))*R*Ta*lift;
    end
    idx = LM(:,e)';
    F(idx) = F(idx)+ F_e;
end
            
        
    
end