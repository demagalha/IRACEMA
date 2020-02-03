% Material Tensor array for usage in Sum Factorization
%% INPUTS
% D - Hooke's Material Law Matrix. (isotropic/anisotropic/etc)
% J - Inverse of Jacobian Matrix @ quadrature points. size = (3,3,m)
% Jmod - Jacobian determinant @ quadrature points.
% i,j,k,l - Integers representing the components to be evaluated
%% OUTPUTS
% Dhat - Column vector representing the entry D_ijkl(x)
function Dhat = MaterialTensor(D,J,Jmod,i,j,k,l) 
    E1 = [1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 0 1 0;];
    E2 = [ 0 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 0; 1 0 0];
    E3 = [ 0 0 0; 0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 0];

    [sz1, sz2] = size(E1);
    E = zeros(sz1,sz2,3);
    E(:,:,1) = E1;
    E(:,:,2) = E2;
    E(:,:,3) = E3;

    [~,~,m] = size(J);
    % Loop over quadrature points
    for x=1:m
        if (cond(J(:,:,x)) > 1e3);
            DF = zeros(3);
        else
            DF = inv(J(:,:,x));
        end
        Ehat1 = E(:,:,i)*DF(k,:)';
        Ehat2 = E(:,:,j)*DF(l,:)';
        % Hooke's Law matrix in parameter domain
        tmp = (Ehat1'*D)*Ehat2;
         Dhat(x) = tmp*Jmod(x);
    end
end