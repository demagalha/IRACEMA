%% Stage one of Sum Factorization and Row Loop assembly
%% INPUTS
% C - 3-dimensional Tensor made by calling MaterialTensor function on the
% three directions of quadrature points. size(C) = [q1,q2,q3,1,1]
% trialfuns - Vector of trial basis functions evaluated at quadrature point k in
% coordinate direction 3. size(trialfuns) = [1,q3]
% testfuns - Vector of test basis functions evaluated ar quadrature point k in
% coordinate direction 3 size(testfuns) = [n3,q3]
%% OUTPUTS
% A - semi-sum-factorized matrix over third dimension. size(A) =
% [n3,m3,q1,q2]
function A = SumFactorizationStageOne(C, trialfuns, testfuns)
    [q1, q2, q3, ~, ~] = size(C);
    [n3, ~] = size(testfuns);
    [m3, ~] = size(trialfuns);
    A = zeros(n3,m3,q1,q2);
    % Looping quadrature points in 3rd direction
    for k=1:q3
        B = testfuns(:,k)*(trialfuns(:,k)');
        % Adding contribution from k to A
        for j = 1:q2 % Loop over q2 quad points
            for i=1:q1 % Loop over q1 quad points
                A(:,:,i,j) = A(:,:,i,j) + C(i,j,k,1,1)*B;
            end
        end
    end        
end