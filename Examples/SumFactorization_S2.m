function A = SumFactorization_S2(C,testfuns,trialfuns)
[~,~,q1,q2] = size(C);
% Loop over second array dimension
for j=1:q2
    B = testfuns(:,j)*(trialfuns(:,j)');
    % Loop over first array dimension
    for i=1:q1
        A(:,:,i) = A(:,:,i) +kron(C(:,:,i,j),B);
    end
end
end