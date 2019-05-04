function Nip = OneBasisFun(p,m,U,i,u)

i = i+1;

if ((i == 0 && u == U(1)) || (i == m-p-1 && u == U(m+1)))
    Nip = 1;
    return
end

if (u < U(i) || u >= U(i+p+1))
    Nip = 0;
    return
end

for j=0:p
    if (u >= U(i+j) && u < U(i+j+1))
        N(j+1) = 1;
    else
        N(j+1) = 0;
    end
end

for k=1:p
    if N(1) == 0
        saved = 0;
    else
        saved = ((u - U(i))*N(1))/(U(i+k)-U(i));
    end
    
    for j=0:(p-k)
        Uleft = U(i+j+1);
        Uright = U(i+j+k+1);
        
        if N(j+2) == 0
            N(j+1) = saved;
            saved = 0;
        else
            temp = N(j+2)/(Uright-Uleft);
            N(j+1) = saved + (Uright - u)*temp;
            saved = (u-Uleft)*temp;
        end
    end
end

Nip = N(1);
end

%U = [0,0,0,1,2,3,4,4,5,5,5]
%u = 5/2
%i = 4
%m = numel(U) or length(U) //numel(U) -1
%p = 2
%ans = 0.1250
%i = 3
%ans = 0.75
            
    

