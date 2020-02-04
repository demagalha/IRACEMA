function Dhat = MaterialTensor(D,X,J,Jmod,i,j,k,l)
E1 = [1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 0 1 0;];
E2 = [ 0 0 0; 0 1 0; 0 0 0; 0 0 1; 0 0 0; 1 0 0];
E3 = [ 0 0 0; 0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 0 0];

[sz1, sz2] = size(E1);
E = zeros(sz1,sz2,3);
E(:,:,1) = E1;
E(:,:,2) = E2;
E(:,:,3) = E3;

[~,~,m] = size(J);
 % Loop over quadrature points in x direction
 for x=1:m
     Ehat1 = E(:,:,i)*J(k,:,x);
     Ehat2 = E(:,:,j)*J(l,:,x);
     % Hooke's Law matrix in parameter domain
     tmp = (Ehat1'*D(X(:,x))*Ehat2);
     Dhat(x) = tmp(1)*Jmod(x);
 end
end
 




