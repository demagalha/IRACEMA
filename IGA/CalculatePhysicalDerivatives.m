function [R, dRdx, d2Rdx2, dudx, J, H] = CalculatePhysicalDerivatives(Model, CollocationPoint)

    u = CollocationPoint(1);
    U = Model.U;
    pu = Model.pu;
    nu = length(U)-pu-1-1;
    su = FindSpanLinear(nu,pu,u,U);
    NU = DersBasisFun(su,u,pu,2,U);
    N = NU(1,:);
    dN = NU(2,:);
    d2N = NU(3,:);

    v = CollocationPoint(2);
    V = Model.V;
    pv = Model.pv;
    nv = length(V)-pv-1-1;
    sv = FindSpanLinear(nv,pv,v,V);
    MV = DersBasisFun(sv,v,pv,2,V);
    M = MV(1,:);
    dM = MV(2,:);
    d2M = MV(3,:);

    w = CollocationPoint(3);
    W = Model.W;
    pw = Model.pw;
    nw = length(W)-pw-1-1;
    sw = FindSpanLinear(nw,pw,w,W);
    LW = DersBasisFun(sw,w,pw,2,W);
    L = LW(1,:);
    dL = LW(2,:);
    d2L = LW(3,:);

    [~, si] = size(NU);
    [~, sj] = size(MV);
    [~, sk] = size(LW);

    clear NU MV LW

    P = Model.get_point_cell;
    weight = Model.weight;
    P = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
    weight = weight(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);

    Q = 0;
    dQdu = 0;
    dQdv = 0;
    dQdw = 0;
    d2Qdu2 = 0;
    d2Qdv2 = 0;
    d2Qdw2 = 0;
    d2Qdudv = 0;
    d2Qdudw = 0;
    d2Qdvdw = 0;
    B = zeros(1,si*sj*sk);
    dBdu = zeros(size(B));
    dBdv = dBdu;
    dBdw = dBdu;
    
    for idx = 1:si*sj*sk
        [i,j,k] = ind2sub([si,sj,sk],idx);
        B(idx) = N(i)*M(j)*L(k);
        dBdu(idx) = dN(i)*M(j)*L(k);
        dBdv(idx) = N(i)*dM(j)*L(k);
        dBdw(idx) = N(i)*M(j)*dL(k);
        Q = Q + B(idx)*weight(idx);
        dQdu = dQdu + dBdu(idx)*weight(idx);
        dQdv = dQdv + dBdv(idx)*weight(idx);
        dQdw = dQdw + dBdw(idx)*weight(idx);

        d2Qdu2 = d2Qdu2 + d2N(i)*M(j)*L(k);
        d2Qdv2 = d2Qdv2 + N(i)*d2M(j)*L(k);
        d2Qdw2 = d2Qdw2 + N(i)*M(j)*d2L(k);
        d2Qdudv = d2Qdudv + dN(i)*dM(j)*L(k);
        d2Qdudw = d2Qdudw + dN(i)*M(j)*dL(k);
        d2Qdvdw = d2Qdvdw + N(i)*dM(j)*dL(k);

    end
    % Basis and parametric derivatives
    R = zeros(1,si*sj*sk);
    d2Rdu2 = zeros(size(R));
    d2Rdv2 = d2Rdu2;
    d2Rdw2 = d2Rdu2;
    d2Rdudv = d2Rdu2;
    d2Rdudw = d2Rdu2;
    d2Rdvdw = d2Rdu2;
    dRdu = zeros(size(R));
    dRdv = zeros(size(R));
    dRdw = zeros(size(R));
    
    for idx=1:si*sj*sk
        [i,j,k] = ind2sub([si,sj,sk],idx);
        R(idx) = B(idx)/Q;
        
        dRdu(idx) = dN(i)*M(j)*L(k) -(R(idx)/Q)*dQdu;
        dRdv(idx) = N(i)*dM(j)*L(k) -(R(idx)/Q)*dQdv;
        dRdw(idx) = N(i)*M(j)*dL(k) -(R(idx)/Q)*dQdw;
     % Parametric second derivatives

        % Main
        d2Rdu2(idx) = d2N(i)*M(j)*L(k) - (1/(Q*Q))*(Q*dRdu(idx) - R(idx)*d2Qdu2);
        d2Rdv2(idx) = N(i)*d2M(j)*L(k) - (1/(Q*Q))*(Q*dRdv(idx) - R(idx)*d2Qdv2);
        d2Rdw2(idx) = N(i)*M(j)*d2L(k) - (1/(Q*Q))*(Q*dRdw(idx) - R(idx)*d2Qdw2);

        % Cross
        d2Rdudv(idx) = dN(i)*dM(j)*L(k) - (1/(Q*Q))*(Q*dRdv(idx) - R(idx)*d2Qdudv);
        d2Rdudw(idx) = dN(i)*M(j)*dL(k) - (1/(Q*Q))*(Q*dRdw(idx) - R(idx)*d2Qdudw);
        d2Rdvdw(idx) = N(i)*dM(j)*dL(k) - (1/(Q*Q))*(Q*dRdw(idx) - R(idx)*d2Qdvdw);

    end

    dR = zeros(3,length(R));
    dR(1,:) = dRdu;
    dR(2,:) = dRdv;
    dR(3,:) = dRdw;

    d2R = zeros(3,3,length(R));
    d2R(1,1,:) = d2Rdu2;
    d2R(1,2,:) = d2Rdudv;
    d2R(1,3,:) = d2Rdudw;
    d2R(2,1,:) = d2Rdudv;
    d2R(2,2,:) = d2Rdv2;
    d2R(2,3,:) = d2Rdvdw;
    d2R(3,1,:) = d2Rdudw;
    d2R(3,2,:) = d2Rdvdw;
    d2R(3,3,:) = d2Rdw2;

    dxdu = zeros(3);
    for idx=1:si*sj*sk
        [i,j,k] = ind2sub([si,sj,sk],idx);
        for xx=1:3
            for yy=1:3
                dxdu(xx,yy) = dxdu(xx,yy) + P{i,j,k}(xx)*dR(yy,idx);
            end
        end
    end
    J = dxdu; % Jacobian
    dUdX = inv(J);
    dRdx = zeros(3,length(R));
    for idx=1:si*sj*sk
        for xx=1:3
            for yy=1:3
                dRdx(xx,idx) = dRdx(xx,idx) + dR(yy,idx)*dUdX(yy,xx);
            end
        end
    end


    dudx = dUdX(1,:);
    dvdx = dUdX(2,:);
    dwdx = dUdX(3,:);
end

