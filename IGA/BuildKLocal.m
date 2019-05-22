function Ke = BuildKLocal(dR_dx,Jmod,D)
nen = length(dR_dx);
B=zeros(6,nen*3);
for i=1:nen*3
B(1,1:3:end) = dR_dx(:,1);
B(2,2:3:end) = dR_dx(:,2);
B(3,3:3:end) = dR_dx(:,3);
B(4,1:3:end) = dR_dx(:,2);
B(4,2:3:end) = dR_dx(:,1);
B(5,1:3:end) = dR_dx(:,3);
B(5,3:3:end) = dR_dx(:,1);
B(6,2:3:end) = dR_dx(:,3);
B(6,3:3:end) = dR_dx(:,2);
Ke = B'*D*B*Jmod;
end