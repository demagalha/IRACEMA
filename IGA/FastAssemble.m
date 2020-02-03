% function [K,M] = FastAssemble(Model, YOUNG, POISSON)
%% Get Weighted Quadrature Rules for U, V and W knot vectors of the model
URules = getWQ(Model.pu,Model.U);
VRules = getWQ(Model.pv,Model.V);
WRules = getWQ(Model.pw,Model.W);
qu = URules.Points;
qv = VRules.Points;
qw = WRules.Points;
%% Get the Jacobian at each quadrature point
[J,Jmod, qpoints] = coordinates_and_jacobian(Model,URules, VRules, WRules);
D = get_matprop_matrix(1, YOUNG, POISSON);
%% Calculate the C_ijkl material law for each quadrature point
for j=1:3
    for i=1:3
        pull_back = @(k,l) MaterialTensor(D,J,Jmod,i,j,k,l);
        K{i,j} = form_local_stiffness(qu,qv,qw,pull_back);
    end
end
% end