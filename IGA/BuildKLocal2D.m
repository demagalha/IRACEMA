function K_local = BuildKLocal2D(dR,Jmod,C)
nen = length(dR);
B = zeros(3,nen*2);
B(1,1:2:end) = dR(:,1);
B(2,2:2:end) = dR(:,2);
B(3,1:2:end) = dR(:,2);
B(3,2:2:end) = dR(:,1);
K_local = B'*C*B*Jmod;
end