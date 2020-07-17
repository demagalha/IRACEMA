function [R, dR, J] = FastBoundaryShape2D(GeometryObj,IntegrationPoint, ...
    global_basis_index,element_local_mapping,element_ranges,element, Boundary)
%   We have to discover first which direction is being integrated.
%   Since the default numbering of Boundaries is Gamma1,2 \in u and
%   Gamma3,4 \in v, we do
    if Boundary > 2
        QUAD_DIRECTION = 1;
    else
        QUAD_DIRECTION = 2;
    end
    FIXED_DIRECTION = setdiff([1 2],QUAD_DIRECTION);
    
    qu = IntegrationPoint(1);
    
    p = GeometryObj.PolynomialOrder;
    pq = p(QUAD_DIRECTION);
    pb = p(FIXED_DIRECTION);
    
    Knots = GeometryObj.KnotVectorCell;
    QK = Knots{QUAD_DIRECTION};
    BK = Knots{FIXED_DIRECTION};
    
    support = global_basis_index(element_local_mapping(:,element),:);
    
    q_range = element_ranges(element,:,QUAD_DIRECTION);
    
    q = ((q_range(2)-q_range(1))*qu +(sum(q_range)))/2; % Parent -> Parametric
    b = (1-mod(Boundary,2)); % Boundary values are either 1 or 0.
    
    sq = FindSpanLinear(length(QK)-pq-1,pq,q,QK);
    sb = FindSpanLinear(length(BK)-pb-1,pb,b,BK);
    
    P = GeometryObj.get_point_cell;
    ind = sub2ind(size(P),support(:,1),support(:,2));
    ActivePoints = P(ind);
    ActivePoints = cell2mat(ActivePoints);
    Weights = ActivePoints(:,4);
    P = ActivePoints(:,1:3);
    
    NQ = DersBasisFun(sq,q,pq,1,QK);
    if b == 1
        NB = DersBasisFun(sb-1,b,pb,1,BK);
    elseif b == 0
        NB = DersBasisFun(sb,b,pb,1,BK);
    end
    
    N = cell(2,1);
    N{QUAD_DIRECTION} = NQ;
    N{FIXED_DIRECTION} = NB;
    
    M = N{2};
    N = N{1};
    B = kron(M(1,:),N(1,:));
    dBdu = kron(M(1,:),N(2,:));
    dBdv = kron(M(2,:),N(1,:));
    
    Q = B*Weights;
    dQdu = dBdu*Weights;
    dQdv = dBdv*Weights;
    
    R = B'.*Weights/Q;
    
    ratios = Weights/(Q*Q);
    
    dRdu = ratios.*(Q*dBdu' -B'*dQdu);
    dRdv = ratios.*(Q*dBdv' -B'*dQdv);
    
    dxdu = sum(P.*dRdu);
    dxdv = sum(P.*dRdv);

    dXdU = [dxdu', dxdv'];
    dUdX = pinv(dXdU);

    dudx = dUdX(:,1)';
    dvdx = dUdX(:,2)';

    dR = dRdu*dudx +dRdv*dvdx;
    tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
    tmp = [squeeze(tmp); 0];
    dQdU = eye(3);
    dQdU(1,1) = tmp(1);
    dQdU(2,2) = tmp(2);
    dQdU(3,3) = tmp(3);
    Jacobian = dXdU(:,1)*dQdU(1,:) + dXdU(:,2)*dQdU(2,:);
    Jacobian = Jacobian(1:3,1:2);
    J = det(Jacobian(1:2,1:2));
end