function [nh,Uh,Qw] = DegreeElevateCurve(n,p,U,Pw,t)

bezalfs = zeros(p+t+1,p+1); %coeficientes de bezier para elevação de grau
bpts(1:p+1) = CPOINT(0,0,0,0,1);          %pontos de controle de bezier de grau p do segmento atual
ebpts(1:p+t+1) = CPOINT(0,0,0,0,1);       %pontos de controle de bezier de grau p+t do segmento atual
Nextbpts(1:p-1) = CPOINT(0,0,0,0,1);
alphas = zeros(p-1);

ss = FindS(U,p);

Qw(1:n + t*(1+ss)+1) = CPOINT(0,0,0,0,1);

m = n+p+1;
ph = p+t;
ph2 = floor(ph/2);

bezalfs(1,1) = 1;
bezalfs(ph+1,p+1) = 1;

for i=1:ph2
    inv = 1/nchoosek(ph,i);
    mpi = min(p,i);
    for j=max(0,i-t):mpi
        bezalfs(i+1,j+1) = inv*nchoosek(p,j)*nchoosek(t,i-j);
    end
end

for i=ph2+1:ph-1
    mpi = min(p,i);
    for j=max(0,i-t):mpi
        bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);
    end
end

mh = ph;
kind = ph+1;
r = -1;
a = p;
b = p+1;
cind = 1;
ua = U(1);
Qw(1) = Pw(1);

for i=0:ph
    Uh(i+1) = ua;
end

for i=0:p
    bpts(i+1) = Pw(i+1);
end

while b < m
    i = b;
    while (b < m && U(b+1) == U(b+2))
        b = b +1;
    end
    mul = b-i+1;
    mh = mh+mul+t;
    ub = U(b+1);
    oldr = r;
    r = p-mul;
    %insert knot U(b+1) r times
    if oldr > 0
        lbz = floor((oldr+2)/2);
    else
        lbz = 1;
    end
    
    if r > 0
        rbz = ph-floor((r+1)/2);
    else
        rbz = ph;
    end
    
    if r > 0
        %Insert knot to get Bezier segment
        numer = ub-ua;
        for k=p:-1:mul+1
            alfs(k-mul) = numer/(U(a+k+1)-ua);
        end
        
        for j=1:r
            save = r-j;
            s = mul+j;
            for k=p:-1:s
                bpts(k+1) = alfs(k-s+1)*bpts(k+1) + (1-alfs(k-s+1))*bpts(k);
            end
            Nextbpts(save+1) = bpts(p+1);
        end
    end
    
    %degree elevate bezier
    for i=lbz:ph
        ebpts(i+1) = CPOINT(0,0,0,0,1);
        mpi = min(p,i);
        for j=max(0,i-t):mpi
            ebpts(i+1) = ebpts(i+1) + bezalfs(i+1,j+1)*bpts(j+1);
        end
    end
    
    if oldr > 1
        first = kind-2;
        last = kind;
        den = ub-ua;
        bet = (ub-Uh(kind))/den;
        for tr=1:oldr-1
            %knot removal loop
            i = first; j = last; kj = j-kind+1;
            while j-i > tr
                if i < cind
                    alf = (ub-Uh(i+1))/(ua-Uh(i+1));
                    Qw(i+1) = alf*Qw(i+1) + (1-alf)*Qw(i);
                end
                if j >= lbz
                    if (j-tr <= kind-ph+oldr)
                        gam = (ub-Uh(j-tr+1))/den;
                        ebpts(kj+1) = gam*ebpts(kj+1) + (1-gam)*ebpts(kj+2);
                    else
                        ebpts(kj+1) = bet*ebpts(kj+1)+(1-bet)*ebpts(kj+2);
                    end
                end
                i = i +1; j = j-1; kj = kj-1;
            end
            first = first-1; last = last+1;
        end
    end
        
        if a ~= p
            for i=0:ph-oldr-1
                Uh(kind+1) = ua; kind = kind + 1;
            end
        end
        
        for j=lbz:rbz
            Qw(cind+1) = ebpts(j+1); cind = cind +1;
        end
        
        if b < m
            for j=0:r-1
                bpts(j+1) = Nextbpts(j+1);
            end
            for j=r:p
                bpts(j+1) = Pw(b-p+j+1);
            end
            a = b; b = b +1; ua = ub;
        else
            
            for i=0:ph
                Uh(kind+i+1) = ub;
            end
        end
end
nh = mh-ph-1;
end
    










