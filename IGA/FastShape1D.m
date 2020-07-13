function [R, dR, J] = FastShape1D(GeometryObject,IntegrationPoint,global_basis_index, element_local_mapping,element)

qu = IntegrationPoint(1);
pu = GeometryObject.pu;
U = GeometryObject.U;

support_start = global_basis_index(element_local_mapping(1,element),:);
ni = support_start(1);
u = ((U(ni+1)+U(ni))*qu +U(ni+1) +U(ni))/2;

pu = GeometryObject.pu;
P = GeometryObject.get_point_cell;
ActivePoints = P(ni-pu:ni);
ActivePoints = reshape(ActivePoints,numel(ActivePoints),1);
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