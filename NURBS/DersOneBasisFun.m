function ders2 = DersOneBasisFun(p,m,U,i,u,n)

i = i+1;

if (u < U(i) || u >= U(i+p+1))
    for k=0:n
        ders(k+1) = 0;  
    end
    ders2 = ders(k+1);
    return
end

for j=0:p
    if (u >= U(i+j) && u < U(i+j+1))
        N(j+1,1) = 1;
    else
        N(j+1,1) = 0;
    end
end

for k=1:p
    if N(1,k) == 0
        saved = 0;
    else
        saved = ((u-U(i))*N(1,k))/(U(i+k)-U(i));
    end
    for j=0:p-k
        Uleft = U(i+j+1);
        Uright = U(i+j+k+1);
        if N(j+2,k) == 0
            N(j+1,k+1) = saved;
            saved = 0;
        else
            tmp = N(j+2,k)/(Uright - Uleft);
            N(j+1,k+1) = saved + (Uright-u)*tmp;
            saved = (u-Uleft)*tmp;
        end
    end
end

ders(1) = N(1,p+1);

for k=1:n
    for j=0:k
        ND(j+1) = N(j+1,p-k+1);
    end
    for jj=1:k
        if ND(1) == 0
            saved = 0;
        else
            saved = ND(1)/(U(i+p-k+jj)-U(i));
        end
        for j=0:k-jj
            Uleft = U(i+j+1);
            Uright = U(i+j+p-k+jj+1);
            if ND(j+2) == 0
                ND(j+1) = (p-k+jj)*saved;
                saved = 0;
            else
                tmp = ND(j+2)/(Uright - Uleft);
                ND(j+1) = (p-k+jj)*(saved-tmp);
                saved = tmp;
            end
        end
    end
    ders(k+1) = ND(1);
end
ders2 = ders(k+1);
end

%p = 2;
%m = numel(U) -1;
%U = [0,0,0,1,2,3,4,4,5,5,5];
%n = 2 <<< a n-ésima derivada, não é igual a n = m-p-1
%u = 5/2
%i = 4
%output 0.5


    
    
    
    
    
    
    
    
    