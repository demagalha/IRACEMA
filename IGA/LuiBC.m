function K1 = LuiBC(SubDomain1,SubDomain2,K1,BoundaryElementData,InvFun2)
[INN, IEN, nel, nen] = SubDomain1.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
for i=1:nel
    LM(:,i) = reshape(ID(:,IEN(:,i)),nen,1);
end
pu1 = SubDomain1.pu;
pv1 = SubDomain1.pv;

U1 = SubDomain1.U;
V1 = SubDomain1.V;

U2 = SubDomain2.U;
V2 = SubDomain2.V;

Knots1 = {U1, V1};
p1 = [pu1, pv1];

% Quadrature Points
   [u, wu] = getGP(p(1));
    [v, wv] = getGP(p(2));
    quadpoints = {u, wu; v, wv};
    N_QUAD_U = length(u);
    N_QUAD_V = length(v);
    N_QUAD = [N_QUAD_U; N_QUAD_V];
    N_ELE_DOF = nen;
elements = BoundaryElementData;
for ee=1:numel(elements)
    e = elements(ee);
    ni = INN(IEN(1,e),:);
    BOUNDARY_DIRECTION = BoundaryData(ee,2);
    QUAD_DIRECTION = setdiff([1,2],parametric_direction);
    q = quadpoints{QUAD_DIRECTION,1};
    w = quadpoints{QUAD_DIRECTION,2};
    Domain = Knots1{QUAD_DIRECTION};
    n = ni(QUAD_DIRECTION);
    h = Domain(n+1) - Domain(n);
    if h < sqrt(eps)
        continue
    end
    BVAL = BoundaryData(ee,3);
    K_e = zeros(prod(p+1));
    for i=1:N_QUAD(QUAD_DIRECTION)
        qq = (h*q(i) +Domain(n+1)+Domain(n))/2;
        ParameterSpace1 = zeros(2,1);
        ParameterSpace1(QUAD_DIRECTION,BOUNDARY_DIRECTION) = [qq, BVAL];
        uu = ParameterSpace1(1);
        vv = ParameterSpace1(2);
        
        % Getting information from the other subdomain
        PhysicalSpace = SubDomain1.eval_point(uu,vv);
        x = PhysicalSpace.x;
        y = PhysicalSpace.y;
        ParameterSpace2 = InvFun2(x,y);
        u2 = ParameterSpace2(1);
        v2 = ParameterSpace2(2);
        
         % Deformation information
        PhysicalSpace = SubDomain2.eval_point(u2,v2);
        z = PhysicalSpace.z;
        
        % Here, we do a numerical derivative.
        hx = sqrt(eps)*x; % Numerical stable x step
        xph = x + hx;
        
        SteppedParameterSpace2 = InvFun2(xph,y);
        uph2 = SteppedParameterSpace2(1);
        vph2 = SteppedParameterSpace2(2);
        SteppedPhysicalSpace = SubDomain2.eval_point(uph2,vph2);
        zph = SteppedPhysicalSpace.z;
        dzdx = (zph - z)/hx;
        
        % Now we repeat for y
        hy = sqrt(eps)*y; % Numerical stable y step
        yph = y + hy;
        SteppedParameterSpace2 = InvFun2(x,yph);
        uph2 = SteppedParameterSpace2(1);
        vph2 = SteppedParameterSpace2(2);
        SteppedPhysicalSpace = SubDomain2.eval_point(uph2,vph2);
        zph = SteppedPhysicalSpace.z;
        dzdy = (zph - z)/hy;
        
        % Gradient Vector 
        gradz = [dzdx; dzdy];
        
        % Normal direction calculation
        n = get_normal_vector(Subdomain1,BOUNDARY_DIRECTION,BVAL);
        n = n/norm(n); % unitary vector
        
        dz = dot(gradz,n);
        if abs(z) > sqrt(eps)
            BETA = -dz/z;
        else
            BETA = 10e4;
        end
        
        BoundaryData = [BOUNDARY_DIRECTION, QUAD_DIRECTION, BVAL];
        [R, dR, J] = BoundaryShape(SubDomain1,q(i),e,IEN,INN,BoundaryData);
        K_e = K_e + abs(J)*w(i)*R*BETA*R';
    end
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + K_e;
end