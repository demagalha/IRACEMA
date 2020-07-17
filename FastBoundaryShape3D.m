function [R, dR, Jmod] = FastBoundaryShape3D(GeometryObj,IntegrationPoint, ...
    global_basis_index,element_local_mapping,element_ranges,element,Boundary)
if Boundary > 4
    QUAD_DIRECTION = [1,2];
elseif (Boundary > 2) && (Boundary < 5)
    QUAD_DIRECTION = [1,3];
else
    QUAD_DIRECTION = [2,3];
end
FIXED_DIRECTION = setdiff([1 2 3], QUAD_DIRECTION);

qu = IntegrationPoint(1);
qv = IntegrationPoint(2);

p = GeometryObj.PolynomialOrder;
pq = p(QUAD_DIRECTION);
pb = p(FIXED_DIRECTION);

Knots = GeometryObj.KnotVectorCell;
QK = Knots(QUAD_DIRECTION);
BK = Knots{FIXED_DIRECTION};

support = global_basis_index(element_local_mapping(:,element),:);

q_range = element_ranges(element,:,QUAD_DIRECTION);
q_range = squeeze(q_range)';

q1 = ((q_range(1,2)-q_range(1,1))*qu +sum(q_range(1,:)))/2; % Parent -> Parametric
q2 = ((q_range(2,2)-q_range(2,1))*qv +sum(q_range(2,:)))/2;
b = 1-mod(Boundary,2);

sq1 = FindSpanLinear(length(QK{1})-pq(1)-1,pq(1),q1,QK{1});
sq2 = FindSpanLinear(length(QK{2})-pq(2)-1,pq(2),q1,QK{2});
sb = FindSpanLinear(length(BK)-pb-1,pb,b,BK);

P = GeometryObj.get_point_cell;
ind = sub2ind(size(P),support(:,1),support(:,2));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

NQ1 = DersBasisFun(sq1,q1,pq(1),1,QK{1});
NQ2 = DersBasisFun(sq2,q2,pq(2),1,QK{2});
NB = DersBasisFun(sb,b,pb,1,BK);

N = cell(3,1);
N(QUAD_DIRECTION) = {NQ1; NQ2};
N{FIXED_DIRECTION} = NB;

L = N{3};
M = N{2};
N = N{1};

B = kron(L(1,:),kron(M(1,:),N(1,:)));
dBdu = kron(L(1,:),kron(M(1,:),N(2,:)));
dBdv = kron(L(1,:),kron(M(2,:),N(1,:)));
dBdw = kron(L(2,:),kron(M(1,:),N(1,:)));

Q = B*Weights;
dQdu = dBdu*Weights;
dQdv = dBdv*Weights;
dQdw = dBdw*Weights;

R = B'.*Weights/Q;
ratios = Weights/(Q*Q);

dRdu = ratios.*(Q*dBdu' -B'*dQdu);
dRdv = ratios.*(Q*dBdv' -B'*dQdv);
dRdw = ratios.*(Q*dBdv' -B'*dQdw);

dxdu = sum(P.*dRdu);
dxdv = sum(P.*dRdv);
dxdw = sum(P.*dRdw);

dXdU = [dxdu', dxdv', dxdw'];
tmp = pinv(dXdU(:,QUAD_DIRECTION));
dUdX = zeros(3);
dUdX(FIXED_DIRECTION,FIXED_DIRECTION) = 1;
dUdX(QUAD_DIRECTION,QUAD_DIRECTION) = tmp(:,QUAD_DIRECTION);

dudx = dUdX(:,1)';
dvdx = dUdX(:,2)';
dwdx = dUdX(:,3)';

dR = dRdu*dudx +dRdv*dvdx +dRdw*dwdx;

tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
tmp = squeeze(tmp)';
dQdU = eye(3);
dQdU(1,1) = tmp(1);
dQdU(2,2) = tmp(2);
dQdU(3,3) = tmp(3);
dQdU(FIXED_DIRECTION,FIXED_DIRECTION) = 1;

Jacobian = dXdU(:,1)*dQdU(1,:) + dXdU(:,2)*dQdU(2,:) +dXdU(:,3)*dQdU(3,:);
Jacobian(FIXED_DIRECTION,FIXED_DIRECTION) = 1;
Jmod = det(Jacobian);
end