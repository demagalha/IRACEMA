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

%wx = [0,0,1,1,2,3,4,5];
%wy = [0,0,1,1,2,3,4,5];
%wz = [0,0,1,1,2,3,4,5];

%p = 2;
%U = [0,0,0,1,2,3,4,4,5,5,5];
%m = numel(U) -1;
%n = m-p-1;
%u = 5/2;
%d = 1;
%CUR = CreateStruct(wx,wy,wz,0,n,p,U)

% CurveDerivsAlg1dot2(CUR.curve.polygon.n,CUR.curve.p,CUR.curve.U,CUR.curve.polygon.P.x,u,d)
% CurveDerivsAlg1dot2(CUR.curve.polygon.n,CUR.curve.p,CUR.curve.U,CUR.curve.polygon.P.y,u,d)
% CurveDerivsAlg1dot2(CUR.curve.polygon.n,CUR.curve.p,CUR.curve.U,CUR.curve.polygon.P.z,u,d)

%retorna uma array CK onde CK(k+1) é a k-ésima derivada
        

