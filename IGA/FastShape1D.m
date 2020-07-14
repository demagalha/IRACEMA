function [R, dR, J] = FastShape1D(GeometryObject,IntegrationPoint, ...
    global_basis_index, element_local_mapping, element_ranges, element)

qu = IntegrationPoint(1);
pu = GeometryObject.pu;
U = GeometryObject.U;

support = global_basis_index(element_local_mapping(:,element),:);

tmp = sum(element_ranges(element,:,:));  % Equivalent to U(ni+1)+U(ni)
u = tmp(1)*(1+qu)/2; % Parent Coordinates -> Parametric Coordinates

su = FindSpanLinear(length(U)-pu-1,pu,u,U);

P = GeometryObject.get_point_cell;
ind = sub2ind(size(P),support(:,1));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

Basis = DersBasisFun(ni-1,u,pu,1,U);
B = Basis(1,:);
dB = Basis(2,:);
Q = B*Weights;
dQ = dB*Weights;

R = B'.*Weights/W;
dRdu = Weights.*(Q*dB'-dQ*B')/(Q*Q);
x = sum((R.*P));

dxdu = P.*dRdu;
J = norm(sum(dxdu))*(U(ni+1)-U(ni))/2;
dR = dRdu/J;
end