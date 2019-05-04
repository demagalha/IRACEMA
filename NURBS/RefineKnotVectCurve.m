function [Ubar, Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)

m = n+p+1;
a = FindSpanLinear(n,p,X(1),U);
b = FindSpanLinear(n,p,X(r+1),U);
b = b+1;

Qw(1:numel(Pw)+numel(X)-1) = CPOINT(0,0,0,1,1);
for j=0:a-p
    Qw(j+1) = Pw(j+1);
end

for j=b-1:n
    Qw(j+r+2) = Pw(j+1);
end

for j=0:a
    Ubar(j+1) = U(j+1);
end

for j=b+p:m
    Ubar(j+r+2) = U(j+1);
end

i = b+p-1;
k = b+p+r;

for j=r:-1:0
    while (X(j+1) <= U(i+1) && i > a)
        Qw(k-p) = Pw(i-p);
        Ubar(k+1) = U(i+1);
        k = k -1;
        i = i -1;
    end
    Qw(k-p) = Qw(k-p+1);
    for L=1:p
        ind = k-p+L;
        alfa = Ubar(k+L+1) - X(j+1);
        if abs(alfa) == 0
            Qw(ind) = Qw(ind+1);
        else
            alfa = alfa/(Ubar(k+L+1)-U(i-p+L+1));
            Qw(ind) = alfa*Qw(ind) + (1-alfa)*Qw(ind+1);
        end    
    end
    Ubar(k+1) = X(j+1);
    k = k-1;
end
end
    