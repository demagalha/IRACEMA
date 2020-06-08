function H1_norm = H1Norm(Model)
L2 = L2Norm(Model);
semi = H1SemiNorm(Model);
H1_norm = sqrt(L2^2 + semi^2);
end