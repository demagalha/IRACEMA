%% Stage two of Sum Factorization and Row Loop assembly
%% INPUTS
% C - Four-dimensional tensor obtained in StageOne. size(C) = [n3,m3,q1,q2]
% trialfuns - Vector of trial basis functions evaluated at quadrature point j in
% coordinate direction 3. size(trialfuns) = [1,q2]
% testfuns - Vector of test basis functions evaluated ar quadrature point j in
% coordinate direction 3 size(testfuns) = [n2,q2]
%% OUTPUTS
% A - semi-sum-factorized matrix over third dimension. size(A) =
% [n3*n2,m3*m2,q1]
function A = SumFactorizationStageTwo(C, trialfuns, testfuns)
    [n3,m3,q1,q2] = size(C);
    [n2,~] = size(testfuns);
    [m2,~]= size(trialfuns);
    A = zeros(n3*n2,m3*m2,q1);
    % Loop over 2nd dimension
    for j=1:q2
        B = testfuns(:,j)*(trialfuns(:,j)');
        % Loop over 1st dimension
        for i=1:q1
            A(:,:,i) = A(:,:,i) +kron(C(:,:,i,j),B);
        end
    end
end