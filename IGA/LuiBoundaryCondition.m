function K =  LuiBoundaryCondition(Omega1,Omega2,K,BoundaryData,InverseFunction,Mode)
[INN, IEN, nel, nen] = Omega1.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
    for i=1:nel
        LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
    end
pu = Omega1.pu;
pv = Omega1.pv;
U = Omega1.U;
V = Omega1.V;
p = [pu; pv];
Knots = {U,V};
% P = Omega1.get_point_cell;
 % Get Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu);
    [v, wv] = getGP(pv);
    quadpoints = {u, wu; v, wv};
    N_QUAD_U = length(u);
    N_QUAD_V = length(v);
    N_QUAD = [N_QUAD_U; N_QUAD_V];
N_ELE_DOF = nen;
% Get the elements of the Boundary
elements = BoundaryData(:,1);
for ee=1:numel(elements)
    e = elements(ee);
    K_e = zeros(N_ELE_DOF);
    M_e = K_e;
    ni = INN(IEN(1,e),:);
    parametric_direction = BoundaryData(ee,2); % Direction of the Boundary
    direction = setdiff([1,2],parametric_direction); % Direction of Integration
    pidx = ni(parametric_direction)-p(parametric_direction):ni(parametric_direction);
    q = quadpoints{direction,1};
    w = quadpoints{direction,2};
    Domain = Knots{direction};
    n = ni(direction);
    h = Domain(n+1) - Domain(n);
    if h < sqrt(eps)
        continue
    end
    bound_val = BoundaryData(ee,3);
    robin = BoundaryData(ee,4);
    K_e = zeros(prod(p+1));
    for i=N_QUAD(direction)
        if parametric_direction == 2
            uu = ((U(n+1)-U(n))*u(i) +U(n+1)+U(n))/2; 
            vv = bound_val;
        elseif parametric_direction == 1
            vv = ((V(n+1)-V(n))*v(i) +V(n+1)+V(n))/2;
            uu = bound_val;
        end
        x = Omega1.eval_point(uu,vv).x;
        y = Omega1.eval_point(uu,vv).y;
        u2 = InverseFunction{1}(x);
        v2 = InverseFunction{2}(y);
        switch Mode
            case 'finite'
                if direction == 1
                    delta = [sqrt(eps)*u2, 0];
                elseif direction == 2
                    delta = [0, sqrt(eps)*v2];
                end
                dx0 = Omega2.eval_point(u2,v2);
                dx1 = Omega2.eval_point(u2+delta(1),v2+delta(2));
                dzdx = (dx1.z - dx0.z)/(dx1.x - dx0.x);
                dzdy = (dx1.z - dx0.z)/(dx1.y - dx0.y);

                % Calculate the normal vector
                BoundaryDirection = [dx1.x-dx0.x; dx1.y-dx0.y];
                normal = [0 -1; 1 0]*(BoundaryDirection/norm(BoundaryDirection));
                Omega2_derivatives = [dzdx, dzdy];
                dz2 = dot(Omega2_derivatives,normal);
                z2 = Omega2.eval_point(u2,v2).z;
                % Computing the "lifting" functions:
                g = z2;
                f = dz2;
                BETA = -f/g;
                if isnan(BETA)
                    BETA = 0;
                end
            case 'cont'
                P2 = Omega2.get_point_cell;
                [INN2, IEN2, ~, ~] = Omega2.get_connectivity;
                qu2 = (2*u2 +U2(nu) +U2(nu+1))/(U2(nu+1)-U(nu));
                qv2 = (2*v2 +V2(nv) +V2(nv+1))/(V2(nv+1)-V(nv));
                [R2, dR2, J2] = Shape2D(Omega,qu2,qv2,e2,P2,IEN2,INN2);
        end
                
                
        qu = q(i);

        Boundary = [parametric_direction, direction, bound_val];
        [R, ~, J] = BoundaryShape(Omega1,qu,e,IEN,INN,Boundary);
%         F_e = F_e + abs(J)*w(i)*R*robin;
        K_e = K_e + abs(J)*w(i)*R*BETA*R';
        %% Robin Boundary Condition
    end
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + K_e;
end
end