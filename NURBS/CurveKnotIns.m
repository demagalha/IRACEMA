function [nq,UQ,Qw] = CurveKnotIns(np,p,UP,Pw,u,k,s,r)

mp = np+p+1;
nq = np+r;

UQ = zeros(1,mp+r+1);

%%%%%%%%%%%%%%carregando o novo knot vector
for i=0:k
    UQ(i+1) = UP(i+1);
end

for i=1:r
    UQ(k+i+1) = u;
end

for i=k+1:mp
    UQ(i+r+1) = UP(i+1);
end

%%%%%%%%%%%%salvando pontos inalterados
for i=0:k-p
    Qw(i+1) = Pw(i+1);
end

for i=k-s:np
    Qw(i+r+1) = Pw(i+1);
end

for i=0:p-s
    Rw(i+1) = Pw(k-p+i+1);
end

%inserindo o nó r vezes
for j=1:r
    L = k-p+j;
    for i=0:p-j-s
        alpha = (u-UP(L+i+1))/(UP(i+k+2)-UP(L+i+1));
        Rw(i+1) = alpha*Rw(i+2) + (1-alpha)*Rw(i+1);
    end
    Qw(L+1) = Rw(1);
    Qw(k+r-j-s+1) = Rw(p-j-s+1);
end

%carregando o resto dos pontos de controle
for i=L+1:k-s
    Qw(i+1) = Rw(i-L+1);
end

end









