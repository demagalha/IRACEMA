function [K, M] = PoissonProblemAssembly(GeometryObj,SOLUTION_DIMENSIONS)
    [global_basis_index, element_local_mapping, element_ranges] = ...
        GetConnectivityArrays(GeometryObj);
    [~, local_matrix] = BuildGlobalLocalMatrices(element_local_mapping, ...
        SOLUTION_DIMENSIONS);
    p = GeometryObj.PolynomialOrder;
    [quad_point_index, weights] = GenerateQuadPoints(p);
    N_QUAD_POINTS = length(weights);
    N_DOF = SOLUTION_DIMENSIONS*max(max(element_local_mapping));
    K = zeros(N_DOF);
    M = K;
    [N_ELE_DOF,N_ELEMENTS] = size(element_local_mapping);
%% Assembly Loops
    for e=1:N_ELEMENTS
        M_e = zeros(N_ELE_DOF);
        K_e = M_e;
        for n=1:N_QUAD_POINTS
            IntegrationPoint = quad_point_index(n,:);
            [R, dR, J] = FastShape(GeometryObj,IntegrationPoint, ...
                global_basis_index, element_local_mapping, element_ranges,e);
                Jmod = abs(J*weights(n));
            K_e = K_e + Jmod*(dR*dR');
            M_e = M_e + Jmod*(R*R');
        end
        idx = local_matrix(:,e)';
        K(idx,idx) = K(idx,idx) + K_e;
        M(idx,idx) = M(idx,idx) + M_e;
    end
end