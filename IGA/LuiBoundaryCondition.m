function [K,M] =  LuiBoundaryCondition(K,M,Boundary,direction,border,Model1,Model2,InverseFunction)
[INN, IEN, nel, nen] = Model1.get_connectivity;

ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
    for i=1:nel
        LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
    end
pu = Model1.pu;
pv = Model1.pv;
U = Model1.U;
V = Model1.V;
P = Model1.get_point_cell;
 % Get Gauss-Legendre Quadrature Points
    [u, wu] = getGP(pu);
    [v, wv] = getGP(pv);
    N_QUAD_U = length(u);
    N_QUAD_V = length(v);
    N_ELE_DOF = nen;

% Find the elements of the Boundary
elements = [];
for i=1:size(IEN,2)
    a = intersect(IEN(:,i),Boundary);
    if isempty(a)
        continue
    else
        elements = [elements; i];
    end
end
if strcmp(direction,'u')
    N_QUAD = N_QUAD_U;
    v = ones(N_QUAD,1);
    wv = v;
    if border == 0
        v = -v;
    end
elseif strcmp(direction,'v')
    N_QUAD = N_QUAD_V;
    u = ones(N_QUAD,1);
    wu = u;
    if border == 0
        u = -u;
    end
else
    error('Invalid direction. Please choose string u or string v as input.')
end

for ee=1:numel(elements)
    e = elements(ee);
    K_e = zeros(N_ELE_DOF);
    M_e = K_e;
    nu = INN(IEN(1,e),1);
    nv = INN(IEN(1,e),2);
    for i=N_QUAD
        [R, dR, J] = Shape2D(Model1,u(i),v(i),e,P,IEN,INN);
        Jmod = abs(J*wu(i)*wv(i));
        %% Robin Boundary Condition
        % We have to find the g_l and f_l to add to our matrices
        % This is done by accessing the virtual boundary on Model2
        % There are several ways to estimate the derivative in Omega2,
        % And we will do a finite difference estimation.
        uu = ((U(nu+1)-U(nu))*u(i) +U(nu+1)+U(nu))/2; 
        vv = ((V(nv+1)-V(nv))*v(i) +V(nv+1)+V(nv))/2;
        x = Model1.eval_point(uu,vv).x;
        y = Model1.eval_point(uu,vv).y;
        u2 = InverseFunction{1}(x);
        v2 = InverseFunction{2}(y);
        perturbation = sqrt(eps);
        if border == 1
            perturbation = -perturbation;
        end
        
        dx0 = Model2.eval_point(u2,v2);
        dx1 = Model2.eval_point(u2,v2+perturbation);
        dzdx = (dx1.z - dx0.z)/(dx1.x - dx0.x);
        
        dy0 = Model2.eval_point(u2,v2);
        dy1 = Model2.eval_point(u2+perturbation,v2);
        dzdy = (dy1.z - dy0.z)/(dy1.y - dy0.y);
        if strcmp(direction,'u')
            dz2 = dzdy;
        elseif strcmp(direction,'v')
            dz2 = dzdx;
        end
        
        z2 = Model2.eval_point(u2,v2).z;
        
        % Computing the "lifting" functions:
        g = z2;
        f = -abs(dz2);
%         K_e = K_e + Jmod*((g)*(dR(:,1)*dR(:,1)' +dR(:,2)*dR(:,2)') -f*(R*R'));
          K_e = K_e -60*Jmod*(f/g)*(R*R');
%         M_e = M_e + Jmod*f/g*(R*R');
    end
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + K_e;
    M(idx,idx) = M(idx,idx) + M_e;
end
end