function [K,F] = LinearElasticityAssemble2D(GeometryObj,YOUNG,POISSON, ...
    LOAD)
SOLUTION_DIMENSIONS = 2;
[global_basis_index, element_local_mapping, element_ranges] = ...
    GetConnectivityArrays(GeometryObj);
[~, local_matrix] = BuildGlobalLocalMatrices(element_local_mapping, ...
    SOLUTION_DIMENSIONS);
p = GeometryObj.PolynomialOrder;
[quad_point_index, weights] = GenerateQuadPoints(p);
N_QUAD_POINTS = length(weights);
N_DOF = SOLUTION_DIMENSIONS*max(max(element_local_mapping));
K = zeros(N_DOF);
F = zeros(N_DOF,1);
[N_ELE_DOF, N_ELEMENTS] = size(element_local_mapping);
N_ELE_DOF = N_ELE_DOF*SOLUTION_DIMENSIONS;

D = (YOUNG/(1-POISSON^2))*[1        POISSON     0; 
                           POISSON      1       0; 
                           0            0   (1-POISSON)/2];
%% Assembly Loops
for e=1:N_ELEMENTS
    K_e = zeros(N_ELE_DOF);
    F_e = zeros(N_ELE_DOF,1);
    for n=1:N_QUAD_POINTS
        IntegrationPoint = quad_point_index(n,:);
        [R, dR, J] = FastShape(GeometryObj,IntegrationPoint, ...
            global_basis_index, element_local_mapping, element_ranges,e);
        Jmod = abs(J*weights(n));

        B = zeros(3,length(dR)*2);
        B(1,1:2:end) = dR(:,1);
        B(2,2:2:end) = dR(:,2);
        B(3,1:2:end) = dR(:,2);
        B(3,2:2:end) = dR(:,1);
        
        K_e = K_e + B'*D*B*Jmod;
        
        x_index = 1:2:length(R)*2;
        y_index = 2:2:length(R)*2;
        F_e(x_index) = F_e(x_index) +R*LOAD(1)*Jmod;
        F_e(y_index) = F_e(y_index) +R*LOAD(2)*Jmod;
    end
    idx = local_matrix(:,e)';
    K(idx,idx) = K(idx,idx) +K_e;
    F(idx) = F(idx) + F_e;
end
        
end