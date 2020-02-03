function A = SumFactorization_S3(C,testfuns,trialfuns)
[~,~,q1] = size(C);
% Loop over first dimension
for i=1:q1
    B = testfuns(:,i)*(trialfuns(:,i)');
    A(:,:) = A(:,:) +kron(C(:,:,i),B);
end
end