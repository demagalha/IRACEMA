function K_local = row_form_stiffness(QuadPointArray, MaterialPullBackFun)
[q1, q2, q3] = size(QuadPointArray);
C = zeros(q1,q2,q3,1,1);
% Loop over partial derivatives of the test fun
for k=1:3
    for l=1:3
        C(:,:,k,l) = MaterialPullBackFun(k,l);
        K_local
end