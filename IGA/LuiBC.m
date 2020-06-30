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
    BOUNDARY_DIRECTION = BoundaryData(ee,2); % Boundary Direction
    QUAD_DIRECTION = setdiff([1,2],parametric_direction); % Integration Direction
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
        PhysicalSpace = SubDomain1.eval_point(uu,vv);
        x = PhysicalSpace.x;
        y = PhysicalSpace.y;
        ParameterSpace2 = InvFun2(x,y);
        
        
end