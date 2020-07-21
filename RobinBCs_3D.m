function [K,F] = RobinBCs_3D(GeometryObj,K,F,BoundaryBasis,Boundary, ...
    boundary_values)
     [global_basis_index, element_local_mapping, element_ranges] = ...
         GetConnectivityArrays(GeometryObj);

directions = [1 2 3];
QUAD_DIRECTION = setdiff(directions,round(Boundary/2));

p = GeometryObj.PolynomialOrder;
Knots = GeometryObj.KnotVectorCell;

pq = p(QUAD_DIRECTION);
Knots = Knots(QUAD_DIRECTION);

[boundary_ranges1, eConn1] = KnotConnectivity(pq1,Knots{1});
[boundary_ranges2, eConn2] = KnotConnectivity(pq2,Knots{2});

boundary_ranges

[NUMBER_OF_ELEMENTS1, ELEMENT_DOFS1] = size(eConn1);
[NUMBER_OF_ELEMENTS2, ELEMENT_DOFS2] = size(eConn2);
NUMBER_OF_ELEMENTS = NUMBER_OF_ELEMENTS1*NUMBER_OF_ELEMENTS2;
ELEMENT_DOFS = ELEMENT_DOFS1*ELEMENT_DOFS2;

[quad_point_index, weights] = GenerateQuadPoints(pq);
N_QUAD_POINTS = length(weights);


Points = GeometryObj.get_point_cell;
Points = Points(:);
Points = cell2mat(Points);
Points = Points(BoundaryBasis,:);

for e=1:NUMBER_OF_ELEMENTS
    F_e = zeros(ELEMENT_DOFS,1);
    K_e = zeros(ELEMENT_DOFS);
    r = boundary_values(e,1);
    beta = boundary_values(e,2);
    b_ranges = boundary_ranges(e,:);
    for i=1:N_QUAD_POINTS
        IntegrationPoint = quad_point_index(n,:);
        [R, ~, J] = FastBoundaryShape(GeometryObj,IntegrationPoint, ...
            QUAD_DIRECTION, b_ranges,Points,e,element_connectivity);
        F_e = F_e + R*r*abs(w(i)*J);
        K_e = K_e +R*beta*R'*abs(w(i)*J);
    end
    idx = BoundaryBasis(eConn(e,:));
    F(idx) = F(idx)+F_e;
    K(idx,idx) = K(idx,idx) + K_e;         
end
end