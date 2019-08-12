function CK = CurveDerivsAlg1dot2(n,p,U,P,u,d)

%calcula a d-ésima derivada de um ponto na curva para u fixo
%placeholder name

du = min(d,p);


for k=p+1:d
    CK(k+1) = 0;
end

    span = FindSpanLinear(n,p,u,U);
    nders = DersBasisFun(span,u,p,du,U);
    
    for k=0:du
        CK(k+1) = 0;
        for j=0:p
            CK(k+1) = CK(k+1) + nders(k+1,j+1)*P(span-p+j+1);
            %P(span-p+j) -> P(j+1)
        end
    end
end

