function C = CCurvePoint2(n,p,U,Pw,u)

span = FindSpanLinear(n,p,u,U);
N = BasisFuns(span,u,p,U);

Cw = CPOINT(0,0,0,0,1);
for j=0:p
    Cw = Cw + N(j+1)*Pw(span-p+j+1);    
end

C = Cw/Cw.w;

end

