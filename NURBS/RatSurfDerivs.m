function SKL = RatSurfDerivs(Aders,wders,d)

for k=0:d
    for L=0:d-k
        v = Aders(k+1,L+1);
        for j=1:L
            v = v - nchoosek(L,j)*wders(1,j+1)*SKL(k+1,L-j+1);
        end
        
        for i=1:k
            v = v - nchoosek(k,i)*wders(i+1,1)*SKL(k-i+1,L+1);
            v2 = 0;
            for j=1:L
                v2 = v2 + nchoosek(L,j)*wders(i+1,j+1)*SKL(k-i+1,L-j+1);
            end
            v = v - nchoosek(k+1,i+1)*v2;
        end
        SKL(k+1,L+1) = v/wders(1,1);
    end
end