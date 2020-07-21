function [K,F] = RobinBCs_2D(GeometryObj,K,F,BoundaryBasis,Boundary, ...
    boundary_values)
     [global_basis_index, element_local_mapping, element_ranges] = ...
         GetConnectivityArrays(GeometryObj);

directions = [1 2];
QUAD_DIRECTION = setdiff(directions,round(Boundary/2));
BOUNDARY_DIRECTION = round(Boundary/2);

p = GeometryObj.PolynomialOrder;
Knots = GeometryObj.KnotVectorCell;

pq = p(QUAD_DIRECTION);
Knot = Knots{QUAD_DIRECTION};
[boundary_ranges, eConn] = KnotConnectivity(pq,Knot);
[NUMBER_OF_ELEMENTS, ELEMENT_DOFS] = size(eConn);

[q, w] = getGP(pq);
Points = GeometryObj.get_point_cell;
Points = Points(:);
Points = cell2mat(Points);
Points = Points(BoundaryBasis,:);

for e=1:NUMBER_OF_ELEMENTS
    F_e = zeros(ELEMENT_DOFS,1);
    K_e = zeros(ELEMENT_DOFS);
    r = boundary_values(e,1);
    beta = boundary_values(e,2);
    B_RANGE = boundary_ranges(e,:);
    for i=1:length(q)
        IntegrationPoint = q(i);
        [R, ~, J] = FastBoundaryShape(GeometryObj,IntegrationPoint, ...
            QUAD_DIRECTION, B_RANGE,Points,e,eConn);
        F_e = F_e + R*r*abs(w(i)*J);
        K_e = K_e +R*beta*R'*abs(w(i)*J);
    end
    idx = BoundaryBasis(eConn(e,:));
    F(idx) = F(idx)+F_e;
    K(idx,idx) = K(idx,idx) + K_e;         
end
end