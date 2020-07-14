function [R, dR, J] = FastShape2D(GeometryObject,IntegrationPoint, ... 
    global_basis_index, element_local_mapping,element_ranges, element)
qu = IntegrationPoint(1);
qv = IntegrationPoint(2);

pu = GeometryObject.pu;
pv = GeometryObject.pv;

U = GeometryObject.U;
V = GeometryObject.V;

support = global_basis_index(element_local_mapping(:,element),:);

tmp = sum(element_ranges(element,:,:)); % Equivalent to U(ni+1)+U(ni)
u = (tmp(1)*(1+qu))/2; % Parent Coordinates -> Parametric Coordinates
v = (tmp(2)*(1+qv))/2;

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

dxdu = P.*dRdu;
dxdv = P.*dRdv;
J = [sum(dxdu); sum(dxdv)];
Jquad = J*J';
dR = dRdu*sqrt(det(inv(Jquad)));
end