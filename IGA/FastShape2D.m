function [R, dR, Jmod] = FastShape2D(GeometryObject,IntegrationPoint, ... 
    global_basis_index, element_local_mapping,element_ranges, element)
qu = IntegrationPoint(1);
qv = IntegrationPoint(2);

pu = GeometryObject.pu;
pv = GeometryObject.pv;

U = GeometryObject.U;
V = GeometryObject.V;

support = global_basis_index(element_local_mapping(:,element),:);

u_range = element_ranges(element,:,1);
v_range = element_ranges(element,:,2);

u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric
v = ((v_range(2)-v_range(1))*qv +(sum(v_range)))/2;

su = FindSpanLinear(length(U)-pu-1,pu,u,U);
sv = FindSpanLinear(length(V)-pv-1,pv,v,V);

P = GeometryObject.get_point_cell;
ind = sub2ind(size(P),support(:,1),support(:,2));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

N = DersBasisFun(su,u,pu,1,U);
M = DersBasisFun(sv,v,pv,1,V);

B = kron(M(1,:),N(1,:));
dBdu = kron(M(1,:), N(2,:));
dBdv = kron(M(2,:), N(1,:));

Q = B*Weights;
dQdu = dBdu*Weights;
dQdv = dBdv*Weights;

R = B'.*Weights/Q;

ratios = Weights/(Q*Q);

dRdu = ratios.*(Q*dBdu' -B'*dQdu);
dRdv = ratios.*(Q*dBdv' -B'*dQdv);

x = sum((R.*P));

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

% Jacobian = [dXdU(1,1)*dQdU(1,1) + dXdU(1,2)*dQdU(2,1), ...
%             dXdU(1,1)*dQdU(1,2) + dXdU(1,2)*dQdU(2,2);
%             dXdU(2,1)*dQdU(1,1) + dXdU(2,2)*dQdU(2,1), ...
%             dXdU(2,1)*dQdU(1,2) + dXdU(2,2)*dQdU(2,2);
%             dXdU(3,1)*dQdU(1,1) + dXdU(3,2)*dQdU(2,1), ...
%             dXdU(3,1)*dQdU(1,2) + dXdU(3,2)*dQdU(2,2)     
%             ];
Jacobian = dXdU(:,1)*dQdU(1,:) + dXdU(:,2)*dQdU(2,:);
Jacobian = Jacobian(1:3,1:2);
Jmod = det(Jacobian(1:2,1:2));

end