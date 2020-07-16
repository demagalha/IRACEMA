function [K, F] = ApplyRobinBCs(GeometryObj,K,F,BoundaryElements,Boundary, ...
    boundary_values)
    [global_basis_index, element_local_mapping, element_ranges] = ...
        GetConnectivityArrays(GeometryObj);
    [~, local_matrix] = BuildGlobalLocalMatrices(element_local_mapping, ...
        SOLUTION_DIMENSIONS);
    p = GeometryObj.PolynomialOrder;
    p(Boundary) = [];
    [quad_point_index, weights] = GenerateQuadPoints(p);

    N_QUAD_POINTS = length(weights);
    [N_ELE_DOF, ~] = size(element_local_mapping);
    N_ELEMENTS = length(BoundaryElements);
    for ee = 1:N_ELEMENTS
        e = BoundaryElements(ee);
        K_e = zeros(N_ELE_DOF);
        F_e = zeros(N_ELE_DOF,1);
        for n=1:N_QUAD_POINTS
            IntegrationPoint = quad_points_index(n,:);
            r = boundary_values(e,1);
            beta = boundary_values(e,2);
            [R, ~, J] = FastBoundaryShape(GeometryObj,IntegrationPoint,...
                global_basis_index,element_local_mapping,element_ranges, ...
                e, Boundary);
            Jmod = abs(J*weights(n));
            F_e = F_e + Jmod*R*r;
            K_e = K_e +Jmod*(R*beta*R');
        end
        idx = local_matrix(:,e)';
        F(idx) = F(idx) + F_e;
        K(idx,idx) = K(idx,idx) +K_e;
    end
end