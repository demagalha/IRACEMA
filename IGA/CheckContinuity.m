%% CheckContinuity Function:
% Checks the continuity of the derivative, via finite differences
% Of two overlapping domains of an eigenvalue Schwarz problem.
%% INPUTS
% Omega1 - Geometry class object which has the boundary
% Omega2 - Geometry class object which is the domain that contains the
% boundary
% BoundaryData - Column array with the elements of the Omega1 Boundary
% InvFunction - InverseFunction to carry information between domains
%% OUTPUTS
% DerivativeData - Array that contains derivative data from the interface
% DerivativeData(1) - Element number on Omega1
% DerivativeData(2) - Element's quadrature point number in linear indexing
% DerivativeData(3) - Quad Point 1st parametric direction in Omega1
% DerivativeData(4) - Quad Point 1st parametric direction in Omega2
% DerivativeData(5) - Quad Point 2nd parametric direction in Omega1
% DerivativeData(6) - Quad Point 2nd parametric direction in Omega2
% DerivativeData(7) - dz/dx in Omega1
% DerivativeData(8) - dz/dx in Omega2
% DerivativeData(9) - dz/dy in Omega1
% DerivativeData(10) - dz/dy in Omega2
%% Function
function DerivativeData = CheckContinuity(Omega1,Omega2,Boundary,InvFunction)
% So, we have to evaluate the derivative continuity
% This is a geometric property of both domains, so we have to check the
% continuity through a given boundary
% This means we have to create a set of points in which we evaluate the
% direcctional derivative. This also allows us to evaluate an error metric
% Naturally, the first choice of points would be the points of quadrature,
% so an element loop, like the one made in applying the boundary conditions
% is to be made.
% We carry this information using the inverse function.
[INN, IEN, nel, nen] = Omega1.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
pu = Omega1.pu;
pv = Omega1.pv;
U = Omega1.U;
V = Omega1.V;
Knots = {U; V};
P = Omega1.get_point_cell;
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

r = 1;
array_size = N_QUAD_U*N_QUAD_V*numel(elements);
DerivativeData = zeros(array_size,6);
for ee=1:numel(elements)
    e = elements(ee);
    nu = INN(IEN(1,e),1);
    nv = INN(IEN(1,e),2);
    count = 1;
    for i=1:N_QUAD_U
        for j=1:N_QUAD_V
        uu = ((U(nu+1)-U(nu))*u(i) +U(nu+1)+U(nu))/2; 
        vv = ((V(nv+1)-V(nv))*v(i) +V(nv+1)+V(nv))/2;
        x = Omega1.eval_point(uu,vv).x;
        y = Omega1.eval_point(uu,vv).y;
        u2 = InvFunction{1}(x);
        v2 = InvFunction{2}(y);
        perturbation = sqrt(eps)*[uu, vv; u2, v2];
        % Calculate Derivatives
        % From Omega 1
        O1_coordinates = Omega1.eval_point(uu,vv);
        h1 = perturbation(1,1);
        h2 = perturbation(1,2);
        O1_h = Omega1.eval_point(uu+h1,vv+h2);
        dz = (O1_h.z - O1_coordinates.z);
        dx = O1_h.x - O1_coordinates.x;
        dy = O1_h.y - O1_coordinates.y;
        O1_dzdx = dz/dx;
        O1_dzdy = dz/dy;
        % From Omega 2
        O2_coordinates = Omega2.eval_point(u2,v2);
        h1 = perturbation(2,1);
        h2 = perturbation(2,2);
        O2_h = Omega2.eval_point(u2+h1,v2+h2);
        dz = (O2_h.z - O2_coordinates.z);
        dx = (O2_h.x - O2_coordinates.x);
        dy = (O2_h.y - O2_coordinates.y);
        O2_dzdx = dz/dx;
        O2_dzdy = dz/dy;
        DerivativeData(r,:) = [x,y,O1_dzdx,O1_dzdy,O2_dzdx,O2_dzdy];
        count = count+1;
        r = r+1;
        end
    end
end