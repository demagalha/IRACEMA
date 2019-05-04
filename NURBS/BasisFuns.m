function N = BasisFuns(i,u,p,U)
 
%cálculo das funções base da B-spline
%entrada i, u, p, U
%saída um vetor N das base
N = zeros(p+1,1); 
N(1) = 1;
 
left = zeros(p+1,1);
right = zeros(p+1,1);

    %como o knot span começa em 0
    %jeito pra burlar isso
    i = i+1;
   
    for j=1:p
    
        left(j+1) = u-U(i+1-j);
        right(j+1) = U(i+j)-u;
        saved = 0;
       
        for r=0:j-1
        
            temp = N(r+1)/(right(r+2) + left(j-r+1));
            N(r+1) = saved + right(r+2)*temp;
            saved = left(j-r+1)*temp;
           
        end
        N(j+1) = saved;
    end
end
 
 
%p = 2;
%U = [0,0,0,1,2,3,4,4,5,5,5];
%u = 5/2;
%i = 4;
%BasisFuns(i,u,p,U)
 
%ans =
%0.12500   0.75000   0.12500
%como no livro