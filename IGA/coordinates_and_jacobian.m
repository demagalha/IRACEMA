function [Basis, DerBasis, qpoints, J, iJ, ControlPoints] = coordinates_and_jacobian(Model,Rules)
    P = Model.get_point_cell;
    U = Model.U;
    pu = Model.pu;
    nu = length(U)-pu-1-1;

    V = Model.V;
    pv = Model.pv;
    nv = length(V)-pv-1-1;
    
    W = Model.W;
    pw = Model.pw;
    nw = length(W)-pw-1-1;
    
    nqu = length(Rules{1}.Points);
    qu = Rules{1}.Points;
    nqv = length(Rules{2}.Points);
    qv = Rules{2}.Points;
    nqw = length(Rules{3}.Points);
    qw = Rules{3}.Points;
    nq = nqu*nqv*nqw;
    qpoints = zeros(nqu*nqv*nqw,3);
    [INN, ~, ~, ~] = Model.get_connectivity;
    for k=1:nqw
        for j=1:nqv
            for i=1:nqu
               idx = sub2ind([nqu, nqv, nqw],i,j,k);
               qpoints(idx,:) = [qu(i), qv(j), qw(k)];
            end
        end
    end
 
    ndofu = length(U)-pu-1;
    ndofv = length(V)-pv-1;
    ndofw = length(W)-pw-1;
    nen = (pu+1)*(pv+1)*(pw+1); % Number of Local Basis Functions
    
    for i=1:nqu*nqv*nqw
        [uu, vv, ww] = ind2sub([nqu,nqv,nqw],i);
        N = Rules{1,1}.Basis(:,uu)';
        dN = Rules{1,1}.DerBasis(:,uu)';
        su = FindSpanLinear(nu,pu,qpoints(i,1),U);
        [~, si] = size(N);
        
        M = Rules{2,1}.Basis(:,vv)';
        dM = Rules{2,1}.DerBasis(:,vv)';
        sv = FindSpanLinear(nv,pv,qpoints(i,2),V); 
        [~, sj] = size(M);
        
        L = Rules{3,1}.Basis(:,ww)';
        dL = Rules{3,1}.DerBasis(:,ww)';
        sw = FindSpanLinear(nw,pw,qpoints(i,3),W);        
        [~, sk] = size(L);
        
        Pq = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
        
        Q = 0;
        dQdu = 0;
        dQdv = 0;
        dQdw = 0;
        B = zeros(si*sj*sk,1);
        dBdu = B;
        dBdv = B;
        dBdw = B;
        for idx=1:si*sj*sk
            [ii, jj, kk] = ind2sub([si,sj,sk],idx);
            B(idx) = N(ii)*M(jj)*L(kk);
            dBdu(idx) = dN(ii)*M(jj)*L(kk);
            dBdv(idx) = N(ii)*dM(jj)*L(kk);
            dBdw(idx) = N(ii)*M(jj)*dL(kk);
            Q = Q + B(idx);
            dQdu = dQdu + dBdu(idx)*Pq{ii,jj,kk}(4);
            dQdv = dQdv + dBdv(idx)*Pq{ii,jj,kk}(4);
            dQdw = dQdw + dBdw(idx)*Pq{ii,jj,kk}(4);
        end
        R = zeros(si*sj*sk,1);
        for idx=1:si*sj*sk
            [ii, jj, kk] = ind2sub([si,sj,sk],idx);
            R(idx) = Pq{ii,jj,kk}(4)*B(idx)/Q;
            dRdu(idx,1) = Pq{ii,jj,kk}(4)/(Q*Q)*(Q*dBdu(idx) - B(idx)*dQdu);
            dRdu(idx,2) = Pq{ii,jj,kk}(4)/(Q*Q)*(Q*dBdv(idx) - B(idx)*dQdv);
            dRdu(idx,3) = Pq{ii,jj,kk}(4)/(Q*Q)*(Q*dBdw(idx) - B(idx)*dQdw);
        end
        
        dxdu = zeros(3);
        for idx=1:si*sj*sk
            [ii, jj, kk] = ind2sub([si,sj,sk],idx);
            for xx=1:3
                for yy=1:3
                    dxdu(xx,yy) = dxdu(xx,yy) +Pq{ii, jj, kk}(xx)*dRdu(idx,yy);
                end
            end
        end
        Jacobian = dxdu;
        J{i} = Jacobian;
        iJ{i} = inv(Jacobian);
        Basis{i} = R;
        DerBasis{i} = dRdu;
        ControlPoints{i} = Pq;
    end
   
end

