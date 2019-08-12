function CK = RatCurveDerivs2(Aders,wders,d)

for k=0:d
    
    v = Aders(k+1);
    for i=1:k
        v = v - nchoosek(k,i)*wders(i+1)*CK(k-i+1);

    end
    CK(k+1) = v/wders(1);
end

end