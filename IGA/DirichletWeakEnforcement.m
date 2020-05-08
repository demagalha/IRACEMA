function F = DirichletWeakEnforcement(Model,BoundaryElements)
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
F = zeros(numel(INN(:,1)),1);
N_ELE_DOF = nen;
for ee=1:length(BoundaryElements)
    e = BoundaryElements(ee,1);
    ni = INN(IEN(1,e),:);
    
   parametric_direction = BoundaryElements(ee,2); % Direction of the Boundary
   pidx = ni(parametric_direction)-p(parametric_direction):ni(parametric_direction);
    direction = mod(parametric_direction,2) +1; % Direction of Integration
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
    for i=1:N_QUAD(direction)
        u = h*q(i) +Domain(n+1)+Domain(n);
        u = u/2;
        N = DersBasisFun(n-1,u,pD,1,Domain);
        if direction == 1
            ActivePoints = P(n-pD:n,pidx);
            aux_fun = zeros(1,p(2)+1);
        elseif direction == 2
            ActivePoints = P(pidx,n-pD:n);
            aux_fun = zeros(1,p(1)+1);
        else
            error('Invalid Parametric Direction')
        end
        ActivePoints = reshape(ActivePoints,numel(ActivePoints),1);
        ActivePoints = cell2mat(ActivePoints);
        if bound_val == 0
            aux_fun(1) = 1;
        elseif bound_val == 1
            aux_fun(end) = 1;
        else
            error('Not a Boundary. Please choose 0 or 1');
        end
        if direction == 1
            Basis = kron(aux_fun,N(1,:));
            derB = kron(aux_fun,N(2,:));
        elseif direction == 2
            Basis = kron(N(1,:),aux_fun);
            derB = kron(N(2,:),aux_fun);
        end
        
        Weights = ActivePoints(:,4);
        Points = ActivePoints(:,1:3);
        
        W = Basis*Weights; % Sum of Basis * Weights
        R = (Basis')./Weights; % Rational Basis Funs
        x = Basis*(Points.*Weights); %x y z coordinates
        
        Jacobian = derB*(Points.*Weights);
        J = norm(Jacobian)*0.5*(h);
        dW = derB*Weights; % Sum of Derivatives * Weights
        dRdu = (1/(W^2))*(W*derB -dW*Basis);
        dR = dRdu/J;
        
        % Computation of penalty terms
        % See Bazilevs, Hughes, 2007 for more details
        gamma = -1; % we choose gamma = +1 for stabilization
        C = 2*h*norm(dRdu)/norm(R); % estimate of C. Can be a bit lower
        % To enforce the inequality, we tweak C a little bit
%         C = 0;

        F_e = F_e + J*w(i)*(gamma*(dR')*lift + (C/h)*R*lift);
    end
    idx = LM(:,e)';
    F(idx) = F(idx) + F_e;
end
            
        
    
end