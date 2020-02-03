function A = SumFactorization_S1(C,testfuns,trialfuns)
% Loop over 3rd dimension
[q1, q2, q3, ~, ~] = size(C);
% Loop over third dimension
for k=1:q3
    B = testfuns(:,k)*(trialfuns(:,k)');
    % Loop over second dimension
    for j=1:q2
        % Loop over first dimension
        for i=1:q1
            A(:,:,i,j) = A(:,:,i,j) +C(i,j,k,1,1)*B;
        end
    end
end