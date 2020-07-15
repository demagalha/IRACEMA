function [R, dR, J] = FastShape1D(GeometryObject,IntegrationPoint, ...
    global_basis_index, element_local_mapping, element_ranges, element)

qu = IntegrationPoint(1);
pu = GeometryObject.pu;
U = GeometryObject.U;

support = global_basis_index(element_local_mapping(:,element),:);

u_range = element_ranges(element,:,1);
u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric

su = FindSpanLinear(length(U)-pu-1,pu,u,U);

P = GeometryObject.get_point_cell;
ind = sub2ind(size(P),support(:,1));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

Basis = DersBasisFun(su,u,pu,1,U);
B = Basis(1,:);
dB = Basis(2,:);
Q = B*Weights;
dQ = dB*Weights;

R = B'.*Weights/W;
dRdu = Weights.*(Q*dB'-dQ*B')/(Q*Q);
% x = sum((R.*P));
dxdu = P.*dRdu;
tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
J = norm(sum(dxdu))*tmp(1);
dR = dRdu/J;
end