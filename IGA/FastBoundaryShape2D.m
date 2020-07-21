function [R, dR, J] = FastBoundaryShape2D(GeometryObj,IntegrationPoint, ...
            QUAD_DIRECTION, B_RANGE,Points,element,element_connectivity)
%   We have to discover first which direction is being integrated.
%   Since the default numbering of Boundaries is Gamma1,2 \in u and
%   Gamma3,4 \in v, we do
    qu = IntegrationPoint(1);
    p = GeometryObj.PolynomialOrder;
    pq = p(QUAD_DIRECTION);
    Knots = GeometryObj.KnotVectorCell;
    QK = Knots{QUAD_DIRECTION};
    q = ((B_RANGE(2)-B_RANGE(1))*qu +(sum(B_RANGE)))/2; % Parent -> Parametric
    sq = FindSpanLinear(length(QK)-pq-1,pq,q,QK);
   
    ind = element_connectivity(element,:);
    ActivePoints = Points(ind,:);
    Weights = ActivePoints(:,4);
    P = ActivePoints(:,1:3);
    
    NQ = DersBasisFun(sq,q,pq,1,QK);
    B = NQ(1,:);
    dBdu = NQ(2,:);
    Q = B*Weights;
    dQdu = dBdu*Weights;
    R = B'.*Weights/Q;
    ratios = Weights/(Q*Q);
    dRdu = ratios.*(Q*dBdu' -B'*dQdu);
      
    dxdu = P.*dRdu;
    tmp = B_RANGE(2)- B_RANGE(1);
    J = norm(sum(dxdu)*tmp(1));
    dR = dRdu/J;
end