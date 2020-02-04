function [K, M, ID] = MultiPatchAssemble(Model,MatPropMatrix,RHO,PatchConnectivity)  
    
    [sz1, ~] = size(PatchConnectivity);
    MS = cell(sz1,1);
    overlapsize = zeros(sz1,1);
    INN = cell(size(Model));
    IEN = INN;
    nel = zeros(size(Model));
    nen = nel;

    % Get Master-Slave data from patches
    for i=1:sz1
        master = PatchConnectivity(i,1);
        slave = PatchConnectivity(i,2);
        mface = PatchConnectivity(i,3);
        sface = PatchConnectivity(i,4);
        MS{i} = MasterSlave(Model{master}, Model{slave}, mface, sface);
        overlapsize(i) = length(MS{i});
    end
 
    % Get INN and IEN vectors from Models
    ID = cell(size(Model));
    for i=1:numel(Model)
        [INN{i}, IEN{i}, nel(i), nen(i)] = Model{i}.get_connectivity;
        ID{i} = reshape(1:max(max(IEN{i}))*3,3,max(max(IEN{i})));
    end
    
    for i=1:sz1
        master = PatchConnectivity(i,1);
        slave = PatchConnectivity(i,2);
       % Update slave's ID vector with correct indexes
         ID{slave}(:,MS{i}(:,2)) = ID{master}(:,MS{i}(:,1));
         maxadd = 1;
        for j=1:length(ID{slave})
           % Check if element is a slave element
           if ismember(j,MS{i}(:,2))
               continue
           elseif maxadd == 1
               % If it is not a slave, start the indexes from the last
               % ID's highest DOF.
               ID{slave}(:,j) = ID{slave-1}(:,end) +3;
               maxadd = 0;
               continue
           end
           % Check if previous element is a slave element
           if ismember(j-1,MS{i}(:,2))
              ID{slave}(:,j) = ID{slave}(:,j-2) +3;
           else
              ID{slave}(:,j) = ID{slave}(:,j-1) +3;
           end
       end
    end

    
    % Change the DOF from slave faces to their correct global number
  
    LM = cell(size(Model));
    for i=1:numel(Model)
        LM{i} = zeros(3*nen(i),nel(i));
        for j=1:nel(i)
            LM{i}(:,j)=reshape(ID{i}(:,IEN{i}(:,j)),3*nen(i),1);
        end
    end
    
   
%     [u{1}, wu{1}] = getGP(Model{1}.pu);
%     [v{1}, wv{1}] = getGP(Model{1}.pv);
%     [w{1}, wew{1}] = getGP(Model{1}.pw);
%     [u{2}, wu{2}] = getGP(Model{2}.pu);
%     [v{2}, wv{2}] = getGP(Model{2}.pv);
%     [w{2}, wew{2}] = getGP(Model{2}.pw);
%     [uu{1}, wuu{1}] = getCG(Model{1}.pu);
%     [vv{1}, wvv{1}] = getCG(Model{1}.pv);
%     [ww{1}, www{1}] = getCG(Model{1}.pw);
%     [uu{2}, wuu{2}] = getCG(Model{2}.pu);
%     [vv{2}, wvv{2}] = getCG(Model{2}.pv);
%     [ww{2}, www{2}] = getCG(Model{2}.pw);

%     N_QUAD_U = length(u);
%     N_QUAD_V = length(v);
%     N_QUAD_W = length(w);

    % Alocate memory for Stiffness and Mass arrays
    K = zeros(max(max(ID{end})));
    M = K;
%     check1 = min(INN);
%     check2 = max(INN);
    D = MatPropMatrix;
    
    %% Assembly
for pp=1:numel(Model)
    pu = Model{pp}.pu;
    pv = Model{pp}.pv;
    pw = Model{pp}.pw;
    U = Model{pp}.U;
    V = Model{pp}.V;
    W = Model{pp}.W;
    P = Model{pp}.get_point_cell;
    N_ELE_DOF = nen(pp)*3;
     % Get Gauss-Legendre Quadrature Points
    [u{pp}, wu{pp}] = getGP(pu);
    [v{pp}, wv{pp}] = getGP(pv);
    [w{pp}, wew{pp}] = getGP(pw);
    % Get Cauchy-Galerkin Quadrature Points
    [uu{pp}, wuu{pp}] = getCG(pu);
    [vv{pp}, wvv{pp}] = getCG(pv);
    [ww{pp}, www{pp}] = getCG(pw);
    % Element on boundary condition
    check1 = min(INN{pp});
    check2 = max(INN{pp});
    for e=1:nel(pp) % Loop Through Elements
        ni = INN{pp}(IEN{pp}(1,e),1);
        nj = INN{pp}(IEN{pp}(1,e),2);
        nk = INN{pp}(IEN{pp}(1,e),3);
        % Check if element has zero measure
        if (U(ni+1) == U(ni) || (V(nj+1) == V(nj)) || (W(nk+1) == W(nk)))
            continue
        end
        K_e = zeros(3*nen(pp),3*nen(pp));
        M_e = K_e;
        cond1 = (ni == check1(1) || nj == check1(2) || nk == check1(3));
        cond2 = (ni == check2(1) || nj == check2(2) || nk == check2(3));
        if (cond1 || cond2)
% % %             If on the boundary, use Gauss-Legendre Quadrature
            qu = u{pp}; qwu = wu{pp}; qv = v{pp}; qwv = wv{pp}; qw = w{pp}; qww = wew{pp};
        else
% %             % If on the interior, use Cauchy-Galerkin Points for Quadrature
            qu = uu{pp}; qwu = wuu{pp}; qv = vv{pp}; qwv = wvv{pp}; qw = ww{pp}; qww = www{pp};
        end
        N_QUAD_U = length(qu);
        N_QUAD_V = length(qv);
        N_QUAD_W = length(qw);
%         F_e = zeros(3*nen,1); Currently disabled.
        for i=1:N_QUAD_U % Loop through U quadrature points
            for j=1:N_QUAD_V % Loop through V quadrature points
                for k=1:N_QUAD_W % Loop through W quadrature points
                    [R, dR, J] = Shape3D(Model{pp},qu(i),qv(j),qw(k),e,P,IEN{pp},INN{pp});
                    Jmod = abs(J*qwu(i)*qwv(j)*qww(k));
                    K_e = K_e + BuildKLocal(dR,Jmod,D);
                    M_e = M_e + BuildMLocal(R,Jmod,RHO);
                end
            end
        end
        % ASSEMBLY
        idx = LM{pp}(:,e);
        for i=1:N_ELE_DOF
            ii = idx(i);
            for j=1:N_ELE_DOF
                jj = idx(j);
                K(ii,jj) = K(ii,jj)+K_e(i,j);
                M(ii,jj) = M(ii,jj)+M_e(i,j);
            end
        end
    end
end
    K = sparse(K);
    M = sparse(M);
end