clearvars
close all
clc

%% Schwarz Method test

% Geometry: One membrane, composed of two overlapping square domains
    a = 1;
    b = 1.5;
    overlap = 0.75;
    
    P1 = [0 0 0 1];
    P2 = [0 a 0 1];
    P3 = [b*(1+overlap)/2 0 0 1];
    P4 = [b*(1+overlap)/2 a 0 1];
    c1 = Geometry('curve',1,[0 0 1 1],{P1 P2});
    c2 = Geometry('curve',1,[0 0 1 1],{P3 P4});
    Omega1 = geo_ruled(c1,c2);

    P1 = [b*(1-overlap)/2 0 0 1];
    P2 = [b*(1-overlap)/2 a 0 1];
    P3 = [b 0 0 1];
    P4 = [b a 0 1];
    c1 = Geometry('curve',1,[0 0 1 1],{P1 P2});
    c2 = Geometry('curve',1,[0 0 1 1],{P3 P4});
    Omega2 = geo_ruled(c1,c2);

clearvars c1 c2 P1 P2 P3 P4

%% Refinement and Patch Conform
% k-refinement
    p = 1; % Number of p-refinements
    Omega1.DegreeElevate(p,1);
    Omega1.DegreeElevate(p,2);
    Omega2.DegreeElevate(p,1);
    Omega2.DegreeElevate(p,2);

    
% Add Conforming knots to both patches
    invNURBS1 = @(x) x/(b*(1+overlap)/2);
    conf1 = invNURBS1(b*(1-overlap)/2);
    Omega1.KnotRefine(conf1,2);
    
    invNURBS2 = @(x) invNURBS1(x- b*(1-overlap)/2);
    conf2 = invNURBS2(b*(1+overlap)/2);
    Omega2.KnotRefine(conf2,2);

    r = 20; % Number of h-refinements
    interval = linspace(0,1,r+2);
    interval = interval(2:end-1);
    Omega1.KnotRefine(interval,1);
    Omega1.KnotRefine(interval,2);
    Omega2.KnotRefine(interval,1);
    Omega2.KnotRefine(interval,2);
%% Assembly Omega 1
  [INN1, IEN1, nel1, nen1] = Omega1.get_connectivity;
    ID1 = reshape(1:max(max(IEN1)),1,max(max(IEN1)));
    LM1 = zeros(nen1,nel1);
    for i=1:nel1
        LM1(:,i) = reshape(ID1(:,IEN1(:,i)),nen1,1);
    end
    % Model parameters
    pu1 = Omega1.pu;
    pv1 = Omega1.pv;
    U1 = Omega1.U;
    V1 = Omega1.V;
    P1 = Omega1.get_point_cell;

    % Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu1);
    [v, wv] = getGP(pv1);

    N_QUAD_U = length(u);
    N_QUAD_V = length(v);

    % Alocate Memory for Stiffness and Force Arrays
    N_DOF1 = numel(INN1)/2;
    K1 = zeros(N_DOF1);
    M1 = zeros(N_DOF1);
    N_ELE_DOF1 = nen1;
    for e=1:nel1
        ni = INN1(IEN1(1,e),1);
        nj = INN1(IEN1(1,e),2);
        % Check if element has measure zero
        if (U1(ni+1) == U1(ni) || V1(nj+1) == V1(nj))
            continue
        end
        K_e = zeros(N_ELE_DOF1);
        M_e = zeros(N_ELE_DOF1);
        for i=1:N_QUAD_U
            for j=1:N_QUAD_V
                [R, dR, J] = Shape2D(Omega1,u(i),v(j),e,P1,IEN1,INN1);
                Jmod = J*wu(i)*wv(j);
                K_e = K_e + Jmod*(dR(:,1)*dR(:,1)' +dR(:,2)*dR(:,2)');
                M_e = M_e+ R*R'*Jmod;
            end
        end
        idx = LM1(:,e)';
        K1(idx,idx) = K1(idx,idx) + K_e;
        M1(idx,idx) = M1(idx,idx) + M_e;
    end
    
%% Assembly Omega 2
      [INN2, IEN2, nel2, nen2] = Omega2.get_connectivity;
    ID2 = reshape(1:max(max(IEN2)),1,max(max(IEN2)));
    LM2 = zeros(nen2,nel2);
    for i=1:nel2
        LM2(:,i) = reshape(ID2(:,IEN2(:,i)),nen2,1);
    end
    % Model parameters
    pu2 = Omega2.pu;
    pv2 = Omega2.pv;
    U2 = Omega2.U;
    V2 = Omega2.V;
    P2 = Omega2.get_point_cell;

    % Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu2);
    [v, wv] = getGP(pv2);

    N_QUAD_U = length(u);
    N_QUAD_V = length(v);

    % Alocate Memory for Stiffness and Force Arrays
    N_DOF2 = numel(INN2)/2;
    K2 = zeros(N_DOF2);
    M2 = zeros(N_DOF2);
    N_ELE_DOF2 = nen2;
    for e=1:nel2
        ni = INN2(IEN2(1,e),1);
        nj = INN2(IEN2(1,e),2);
        % Check if element has measure zero
        if (U2(ni+1) == U2(ni) || V2(nj+1) == V2(nj))
            continue
        end
        K_e = zeros(N_ELE_DOF2);
        M_e = zeros(N_ELE_DOF2);
        for i=1:N_QUAD_U
            for j=1:N_QUAD_V
                [R, dR, J] = Shape2D(Omega2,u(i),v(j),e,P2,IEN2,INN2);
                Jmod = J*wu(i)*wv(j);
                K_e = K_e + Jmod*(dR(:,1)*dR(:,1)' +dR(:,2)*dR(:,2)');
                M_e = M_e+ R*R'*Jmod;
            end
        end
        idx = LM2(:,e)';
        K2(idx,idx) = K2(idx,idx) + K_e;
        M2(idx,idx) = M2(idx,idx) + M_e;
    end
%% Boundary Conditions
% Find boundaries
    b1 = find(INN1(:,1) == 1);
    b2 = find(INN1(:,2) == 1);
    b3 = find(INN1(:,1) == max(INN1(:,1)));
    b4 = find(INN1(:,2) == max(INN1(:,2)));
    Gamma1 = union(b1,b2);
    Gamma1 = union(Gamma1,b3);
    % Gamma 1 is composed of the boundaries made by
    % the first and last basis function of the u parametric direction
    % and by the first basis function of the v parametric direction
    
    Gamma_O1 = b4;
    
    % The Gamma OVERLAP 1 is composed by the last basis function of the v
    % paramatric direction
    
    b1 = find(INN2(:,1) == 1);
    b2 = find(INN2(:,2) == 1);
    b3 = find(INN2(:,1) == max(INN2(:,1)));
    b4 = find(INN2(:,2) == max(INN2(:,2)));
    Gamma2 = union(b1,b3);
    Gamma2 = union(Gamma2,b4);
    % Gamma 2 is composed by the boundaries made by the first and last
    % basis functions of the u paramatric domain and by the last basis
    % function of the v parametric direction
    
    Gamma_O2 = b2;
    
%     Gamma OVERLAP 2 is composed by the first basis function of the v
%     parametric direction
% Gamma1 and Gamma2 restrictions
    constNod1 = reshape(ID1(:,Gamma1),numel(ID1(:,Gamma1)),1);
    O1Nod = reshape(ID1(:,Gamma_O1),numel(ID1(:,Gamma_O1)),1);
    constNod2 = reshape(ID2(:,Gamma2),numel(ID2(:,Gamma2)),1);
    O2Nod = reshape(ID2(:,Gamma_O2),numel(ID2(:,Gamma_O2)),1);

%     K1(constNod1,constNod1) = 1e30*eye(length(constNod1));
%     K2(constNod2,constNod2) = 1e30*eye(length(constNod2));
    KK1 = K1;
    KK2 = K2;
    MM1 = M1;
    MM2 = M2;
    KK1(constNod1,:) = [];
    KK1(:,constNod1) = [];
    KK2(constNod2,:) = [];
    KK2(:,constNod2) = [];
    MM1(constNod1,:) = [];
    MM1(:,constNod1) = [];
    MM2(constNod2,:) = [];
    MM2(:,constNod2) = [];
    [aa1,O1] = eigs(KK1,MM1,1,'sm');
    [aa2,O2] = eigs(KK2,MM2,1,'sm');
    if abs(min(aa1)) > abs(max(aa1))
        aa1 = -aa1;
    end
%     aa1 = aa1/max(aa1);
    if abs(min(aa2)) > abs(max(aa2))
        aa2 = -aa2;
    end
%     aa2 = aa2/max(aa2);
    omega1 = sqrt(diag(O1));
    omega2 = sqrt(diag(O2));
    constNod1 = sort(constNod1,'ascend');
    for i=1:numel(constNod1)
        if constNod1(i) == 1
            aa1 = [0; aa1];
        elseif constNod1(i) == length(K1)
            aa1 = [aa1; 0];
        else
            aa1 = [aa1(1:constNod1(i)-1); 0; aa1(constNod1(i):end)];
        end
    end
    constNod2 = sort(constNod2,'ascend');
    for i=1:numel(constNod2)
        if constNod2(i) == 1
            aa2 = [0; aa2];
        elseif constNod2(i) == length(K2)
            aa2 = [aa2; 0];
        else
            aa2 = [aa2(1:constNod2(i)-1); 0; aa2(constNod2(i):end)];
        end
    end     
    Omega1Deformed = VisualizeModes(Omega1,aa1,ID1);
    Omega2Deformed = VisualizeModes(Omega2,aa2,ID2);
    old1 = Omega1Deformed;
    old2 = Omega2Deformed;
%     figure(1)
%     subplot(3,2,1)
%     Omega1.plot_geo;
%     hold on
%     Omega2.plot_geo;
%     title('Configuração Inicial, Overlap 3/4')
%     alpha(0.9)
%     subplot(3,2,2)
%     Omega1Deformed{1}.plot_geo;
%     hold on
%     Omega2Deformed{1}.plot_geo;
%     str = 'Autovalor Inicial'; 
%     title(str);
%     alpha(0.9)
%% Overlapping Schwarz using Robin Conditions (Lui, 2000)
% Find the elements which compose the boundaries
% e_O1 = Elements of Gamma_O1 that are on Omega 2
e_O1 = find(ismember(IEN1,Gamma_O1));
[~, e_O1] = ind2sub(size(IEN1),e_O1);
e_O1 = unique(e_O1);

% e_O2 = Elements of Gamma_O2 that are on Omega 1
e_O2 = find(ismember(IEN2,Gamma_O2));
[~, e_O2] = ind2sub(size(IEN2),e_O2);
e_O2 = unique(e_O2);
                     
    [u1, wu1] = getGP(pu1);
    [v1, wv1] = getGP(pv1);
    [u2, wu2] = getGP(pu2);
    [v2, wv2] = getGP(pv2);

    N_QUAD_U1 = length(u1);
    N_QUAD_V1 = length(v1);
    N_QUAD_U2 = length(u2);
    N_QUAD_V2 = length(v2);
    iter = 0;
    plot = 0;
str1 = '\omega^1_n = ';
str2 = '\omega^2_n = ';
tol = 1;
% iter = 0;
while iter < 4
for ee=1:numel(e_O1)
    e1 = e_O1(ee);
    ni = INN1(IEN1(1,e1),1);
    K1_e = zeros((pu1+1)*(pv1+1));
    M1_e = K1_e;
    % Update Matrix 1 with Robin Conditions
     for i=N_QUAD_U1
         u1 = ((U1(ni+1) - U1(ni))*u(i) +U1(ni+1)+U1(ni))/2;
         Point = Omega1Deformed{1}.eval_point(u1,1);
         xx = Point.x;
         v2 = invNURBS2(xx);
         u2 = u1;
         
         su1 = FindSpanLinear(length(U1)-pu1-1,pu1,u1,U1);
         sv1 = FindSpanLinear(length(V1)-pv1-2,pv1,1,V1);
         
         su2 = FindSpanLinear(length(U2)-pu2-1,pu2,u2,U2);
         sv2 = FindSpanLinear(length(V2)-pv2-1,pv2,v2,V2);

         P1 = Omega1Deformed{1}.get_point_cell;
         NU = DersBasisFun(su1,u1,pu1,1,U1);
         NV = DersBasisFun(sv1,1,pv1,1,V1);
         R1 = kron(NV(1,:),NU(1,:))';
         dR1 = [kron(NV(1,:),NU(2,:))', kron(NV(2,:),NU(1,:))'];
         B1 = P1(su1-pu1+1:su1+1,sv1-pv1+1:sv1+1);
         B1 = reshape(B1,[numel(B1) 1]);
         B1 = cell2mat(B1);
         B1 = B1(:,1:3).*B1(:,4); % Multiplying by weights
         
         R1x = B1.*R1;
         dxdu1 = zeros(2);
         for loc=1:length(R1x)
             for xx=1:2
                 for yy=1:2
                     dxdu1(xx,yy) = dxdu1(xx,yy) +B1(loc,xx)*dR1(loc,yy);
                 end
             end
         end
         du1dx = inv(dxdu1);
         dR1x = zeros(size(dR1));
         for loc=1:length(R1x)
             for xx=1:2
                 for yy=1:2
                     dR1x(loc,xx) = dR1x(loc,xx) + dR1(loc,yy)*du1dx(yy,xx);
                 end
             end
         end
    
         J1 = abs(det(du1dx)*0.5*(U1(ni+1)-U1(ni)));
         % To correctly enforce Robin conditions on Gamma12, we have
         % to evaluate Omega2's u2 and du2 at these points we are providing
         
         P2 = Omega2Deformed{1}.get_point_cell;
         NU = DersBasisFun(su2,u2,pu2,1,U2);
         NV = DersBasisFun(sv2,v2,pv2,1,V2);
         R2 = kron(NV(1,:),NU(1,:))';
         dR2 = [kron(NV(1,:),NU(2,:))', kron(NV(2,:),NU(1,:))'];
         B2 = P2(su2-pu2+1:su2+1,sv2-pv2+1:sv2+1);
         B2 = reshape(B2,[numel(B2) 1]);
         B2 = cell2mat(B2);
         B2 = B2(:,1:3).*B2(:,4); % Multiplying by weights
         
         R2x = B2.*R2;
         
         dxdu2 = zeros(2);
         for loc=1:length(R2x)
             for xx=1:2
                 for yy=1:2
                     dxdu2(xx,yy) = dxdu2(xx,yy) +B2(loc,xx)*dR2(loc,yy);
                 end
             end
         end
         du2dx = inv(dxdu2);
         dR2x = zeros(size(dR2));
         for loc=1:length(R2x)
             for xx=1:2
                 for yy=1:2
                     dR2x(loc,xx) = dR2x(loc,xx) + dR2(loc,yy)*du2dx(yy,xx);
                 end
             end
         end
    
         J2 = abs(det(du2dx));
         % Robin Conditions on Omega1
            % Boundary varies in y direction, so we take the product:
            du1dy = sum(dR1x(:,2));
            du1dx = sum(dR1x(:,1));
            du2dy = sum(dR2x(:,2));
            du2dx = sum(dR2x(:,1));
            beta1 = (sum(R2x(:,3)));
            beta2 = du2dy;
            % r is zero, so we won't compute
            % Dirichlet Conditions on Omega2
%             g =  x2(3);
%             K1_e = K1_e +1e6*(NU(1,:)'*(x-x2))*J1*wu1(i);
            K1_e = K1_e +((beta1-1)*(dR1x(:,1)*dR1x(:,1)' +dR1x(:,2)*dR1x(:,2)') -beta2*(R1*R1'))/J1*wu1(i);
            M1_e = M1_e +((beta1-1)*(R1*R1'));
      end
        index = LM1(:,e1);
        K1(index,index) = K1(index,index) +K1_e; 
        M1(index,index) = M1(index,index) +M1_e;
end
    KK1 = K1;
    MM1 = M1;
    KK1(constNod1,:) = [];
    KK1(:,constNod1) = [];
    MM1(constNod1,:) = [];
    MM1(:,constNod1) = [];

    [aa1,O1] = eigs(KK1,MM1,1,'sm');
    if abs(min(aa1)) > abs(max(aa1))
        aa1 = -aa1;
    end
    omega1 = sqrt(diag(O1));
    constNod1 = sort(constNod1,'ascend');
    for i=1:numel(constNod1)
        if constNod1(i) == 1
            aa1 = [0; aa1];
        elseif constNod1(i) == length(K1)
            aa1 = [aa1; 0];
        else
            aa1 = [aa1(1:constNod1(i)-1); 0; aa1(constNod1(i):end)];
        end
    end
    aa2 = aa2/max(aa2);
    aa1 = aa1/max(aa1);
    Omega1Deformed = VisualizeModes(Omega1Deformed{1},aa1,ID1);
    Omega2Deformed = VisualizeModes(Omega2Deformed{1},aa2,ID2);

    figure(1)
    iter = iter+1;
    subplot(2,2,iter)
    Omega1Deformed{1}.plot_geo('coarse',0,0,[0 1],[0 1]);
    hold on
    Omega2Deformed{1}.plot_geo('coarse',0,0,[0 1],[conf2 1]);
    title(strcat(str2,num2str(omega2)));
    alpha(0.9)
for ee=1:numel(e_O2)
    e2 = e_O2(ee);
    ni = INN2(IEN2(1,e2),1);
    K2_e = zeros((pu2+1)*(pv2+1));
    M2_e = K2_e;
    % Update Matrix 1 with Robin Conditions
     for i=N_QUAD_U2
         u2 = ((U2(ni+1) - U2(ni))*u(i) +U2(ni+1)+U2(ni))/2;
         Point = Omega2Deformed{1}.eval_point(u2,0);
         xx = Point.x;
         v1 = invNURBS2(xx);
         u1 = u2;
         
         su2 = FindSpanLinear(length(U2)-pu2-1,pu2,u2,U2);
         sv2 = FindSpanLinear(length(V2)-pv2-2,pv2,1,V2);
         
         su1 = FindSpanLinear(length(U1)-pu1-1,pu1,u1,U1);
         sv1 = FindSpanLinear(length(V1)-pv1-1,pv1,v1,V1);

         P2 = Omega2Deformed{1}.get_point_cell;
         NU = DersBasisFun(su2,u2,pu2,1,U2);
         NV = DersBasisFun(sv2,1,pv2,1,V2);
         R2 = kron(NV(1,:),NU(1,:))';
         dR2 = [kron(NV(1,:),NU(2,:))', kron(NV(2,:),NU(1,:))'];
         B2 = P2(su2-pu2+1:su2+1,sv2-pv2+1:sv2+1);
         B2 = reshape(B2,[numel(B2) 1]);
         B2 = cell2mat(B2);
         B2 = B2(:,1:3).*B2(:,4); % Multiplying by weights
         
         R2x = B2.*R2;
         dxdu2 = zeros(2);
         for loc=1:length(R2x)
             for xx=1:2
                 for yy=1:2
                     dxdu2(xx,yy) = dxdu2(xx,yy) +B2(loc,xx)*dR2(loc,yy);
                 end
             end
         end
         du2dx = inv(dxdu2);
         dR2x = zeros(size(dR2));
         for loc=1:length(R2x)
             for xx=1:2
                 for yy=1:2
                     dR2x(loc,xx) = dR2x(loc,xx) + dR2(loc,yy)*du2dx(yy,xx);
                 end
             end
         end
    
         J2 = abs(det(du2dx)*0.5*(U2(ni+1)-U2(ni)));
         % To correctly enforce Robin conditions on Gamma12, we have
         % to evaluate Omega2's u2 and du2 at these points we are providing
         
         P1 = Omega1Deformed{1}.get_point_cell;
         NU = DersBasisFun(su1,u1,pu1,1,U1);
         NV = DersBasisFun(sv1,v1,pv1,1,V1);
         R1 = kron(NV(1,:),NU(1,:))';
         dR1 = [kron(NV(1,:),NU(2,:))', kron(NV(2,:),NU(1,:))'];
         B1 = P1(su2-pu2+1:su2+1,sv2-pv2+1:sv2+1);
         B1 = reshape(B1,[numel(B1) 1]);
         B1 = cell2mat(B1);
         B1 = B1(:,1:3).*B1(:,4); % Multiplying by weights
         
         R1x = B1.*R1;
         
         dxdu1 = zeros(2);
         for loc=1:length(R1x)
             for xx=1:2
                 for yy=1:2
                     dxdu1(xx,yy) = dxdu1(xx,yy) +B1(loc,xx)*dR1(loc,yy);
                 end
             end
         end
         du1dx = inv(dxdu1);
         dR1x = zeros(size(dR1));
         for loc=1:length(R1x)
             for xx=1:2
                 for yy=1:2
                     dR1x(loc,xx) = dR1x(loc,xx) + dR1(loc,yy)*du1dx(yy,xx);
                 end
             end
         end
    
         J1 = abs(det(du1dx));
         % Robin Conditions on Omega1
            % Boundary varies in y direction, so we take the product:
            du1dy = sum(dR1x(:,2));
            du1dx = sum(dR1x(:,1));
            du2dy = sum(dR2x(:,2));
            du2dx = sum(dR2x(:,1));
            beta1 = (sum(R1x(:,3)));
            beta2 = -(du1dy);            % r is zero, so we won't compute
            % Dirichlet Conditions on Omega2
%             g =  x2(3);
%             K1_e = K1_e +1e6*(NU(1,:)'*(x-x2))*J1*wu1(i);
            K2_e = K2_e +((beta1-1)*(dR2x(:,1)*dR2x(:,1)' +dR2x(:,2)*dR2x(:,2)') -beta2*(R2*R2'))/J2*wu2(i);
            M2_e = M2_e +(beta1-1)*(R2*R2');
      end
        index = LM2(:,e2)';
        K2(index,index) = K2(index,index) +K2_e;
        M2(index,index) = M2(index,index) +M2_e;
end
    KK2 = K2;
    MM2 = M2;
    KK2(constNod2,:) = [];
    KK2(:,constNod2) = [];
    MM2(constNod2,:) = [];
    MM2(:,constNod2) = [];
    [aa2,O2] = eigs(KK2,MM2,1,'sm');
    if abs(min(aa2)) > abs(max(aa2))
        aa2 = -aa2;
    end
    omega2 = sqrt(diag(O2));
    constNod2 = sort(constNod2,'ascend');
    for i=1:numel(constNod2)
        if constNod2(i) == 1
            aa2 = [0; aa2];
        elseif constNod2(i) == length(K2)
            aa2 = [aa2; 0];
        else
            aa2 = [aa2(1:constNod2(i)-1); 0; aa2(constNod2(i):end)];
        end
    end
    aa1 = aa1/max(aa1);
    aa2 = aa2/max(aa2);
    Omega2Deformed = VisualizeModes(Omega2Deformed{1},aa2,ID2);
    Omega1Deformed = VisualizeModes(Omega1Deformed{1},aa1,ID1);
    iter = iter+1;
    figure(1)
    subplot(2,2,iter)
    Omega1Deformed{1}.plot_geo('coarse',0,0,[0 1],[0 conf1]);
    hold on
    Omega2Deformed{1}.plot_geo('coarse',0,0,[0 1],[0 1]);
    title(strcat(str2,num2str(omega2)));
    alpha(0.9)
tol = abs(omega1-omega2);
end
% figure(2)
% subplot(2,2,1)
% Omega1Deformed{1}.plot_geo('coarse',0,0,[0 1],[0 1]);
% % title(strcat(str1,num2str(omega1)));
% hold on
% subplot(2,2,2)
% Omega2Deformed{1}.plot_geo('coarse',0,0,[0 1],[0 1]);
% alpha([0.9])
% title(strcat(str1,num2str(omega2)));
% sunplot(2,2,3)


%     clearvars -except V1 V2 ID1 ID2 Omega1 Omega2
%     save('SchwarzMembrane.mat')