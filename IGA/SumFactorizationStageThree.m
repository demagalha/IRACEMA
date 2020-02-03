%% Stage three of Sum Factorization and Row Loop assembly
%% INPUTS
% C - Three-dimensional tensor obtained in StageTwo. size(C) = [n3*n2,m3*m2,q1]
% trialfuns - Vector of trial basis functions evaluated at quadrature point j in
% coordinate direction 3. size(trialfuns) = [1,q1]
% testfuns - Vector of test basis functions evaluated ar quadrature point j in
% coordinate direction 3 size(testfuns) = [n1,q1]
%% OUTPUTS
% A - semi-sum-factorized matrix over third dimension. size(A) =
% [n3*n2*n1,m3*m2*m1]
function A = SumFactorizationStageThree(C,trialfuns,testfuns)
    [n32,m32,q1] = size(C);
    [n1, ~] = size(testfuns);
    [m1, ~] = size(trialfuns);
    A = zeros(n32*n1,m32*m1);
    % Loop over first array dimension
    for i=1:q1
        B = testfuns(:,i)*(trialfuns(:,i)');
        A(:,:) = A(:,:) +kron(C(:,:,i),B);
    end
end