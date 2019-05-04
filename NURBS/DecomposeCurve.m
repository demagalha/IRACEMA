function [nb, Qw] = DecomposeCurve(n,p,U,Pw)

m = n+p+1;
a = p;
b = p+1;
nb = 0;

ind =1;
indmax = length(Pw);
while(ind ~= indmax)
    ind = ind+p;
    nb = nb+1;
end
if(mod(n+1,p+1) == 0)
Qw(1:nb,1:p+1) = CPOINT(0,0,0,0,1);
else
    Qw(1:nb+1,1:p+1) = CPOINT(0,0,0,0,1);
end
nb = 0;

for i=0:p
    Qw(nb+1,i+1) = Pw(i+1);
end

while b < m
    i = b;
    while (b < m && U(b+2) == U(b+1))
        b = b + 1;
    end
    mult = b-i+1;
    if mult < p
        numer = U(b+1) - U(a+1);
        for j=p:-1:mult+1
            alphas(j-mult) =  numer/(U(a+j+1)-U(a+1));
        end
            r = p-mult;
            for j=1:r
                save = r-j;
                s = mult+j;
                for k=p:-1:s
                    alpha = alphas(k-s+1);
                    Qw(nb+1,k+1) = alpha*Qw(nb+1,k+1) + (1-alpha)*Qw(nb+1,k);
                end
                if b < m
                    Qw(nb+2,save+1) = Qw(nb+1,p+1);
                end
            end
    end
    nb = nb + 1;
    if b < m
        for i=p-mult:p
            Qw(nb+1,i+1) = Pw(b-p+i+1);
        end
        a = b;
        b = b +1;
    end
end
end