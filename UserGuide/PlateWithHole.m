clearvars
cd ..
run startup.m
cd UserGuide
clc

a = 1;
L = 4;
Omega = BuildHoledPlate;
YOUNG = 1e5;
POISSON = 0.3;
LOAD = [0 0];
Tx = -10;

% Refinement
% pp = 2;
% Omega.DegreeElevate(pp,1);
% Omega.DegreeElevate(pp,2);
h = 2^6;
interval = linspace(0,1,h+2);
interval = setdiff(interval,[0 1]);
Omega.KnotRefine(interval,1);
Omega.KnotRefine(interval,2);


[K,F] = LinearElasticityAssemble2D(Omega,YOUNG,POISSON,LOAD);

Boundaries = GetBoundaries(Omega);

% One can also use Boundaries = GetBoundaries and Basis =  Boundaries{3}
% as the Boundaries where the exact stresses are put. 
% The way this Holed Plate is constructed, the second parametric direction
% at 0 defines the whole left and uppermost boundaries.
% This is because of the way the surface was ruled between the curves
% There are plenty of ways of constructing this Plate. Hughes' approach is
% different than mine.

%% Neumann Boundary Conditions for the Plate with Hole
[global_basis_index, element_local_mapping, element_ranges] = ...
         GetConnectivityArrays(Omega);
[global_id, ~] = BuildGlobalLocalMatrices(element_local_mapping,2);
Boundary = 4;
BoundaryBasis = Boundaries{4};
directions = [1 2];
QUAD_DIRECTION = setdiff(directions,round(Boundary/2));

p = Omega.PolynomialOrder;
Knots = Omega.KnotVectorCell;
pq = p(QUAD_DIRECTION);
Knot = Knots{QUAD_DIRECTION};
[boundary_ranges, eConn] = KnotConnectivity(pq,Knot);
[NUMBER_OF_ELEMENTS, ELEMENT_DOFS] = size(eConn);
[q, w] = getGP(pq);
Points = Omega.get_point_cell;
Points = Points(:);
% Points = cell2mat(Points);
% Points = Points(BoundaryBasis,:);
for e=1:NUMBER_OF_ELEMENTS
    F_e = zeros(ELEMENT_DOFS*2,1);
    B_RANGE = boundary_ranges(e,:);
    for i=1:length(q)
        IntegrationPoint = q(i);
        [R, ~, J] = FastBoundaryShape(Omega,IntegrationPoint, ...
            QUAD_DIRECTION, B_RANGE, Points, e, eConn);
        Jmod = abs(J*w(i));
        % Exact stress at point
        qq = ((B_RANGE(2) - B_RANGE(1))*IntegrationPoint +(sum(B_RANGE)))/2;
        physical_coord = Omega.eval_point(qq,0);
        x = physical_coord.x;
        y = physical_coord.y;
        
        r = sqrt(x^2 +y^2);
        theta = atan(y/x);
        ratio = a/r;
        r2 = ratio^2;
        r4 = ratio^4;
        c2 = cos(2*theta);
        c4 = cos(4*theta);
        s2 = sin(2*theta);
        
        sigma_rr = ((Tx)/2)*(1-r2) +((Tx)/2)*(1-4*r2 +3*r4)*c2;
        sigma_tt = ((Tx)/2)*(1+r2) -((Tx)/2)*(1 +3*r4)*c2;
        sigma_rt = -((Tx)/2)*(1 +2*r2 -3*r4)*s2;
        
        T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        S = T*[sigma_rr sigma_rt; sigma_rt sigma_tt];
        
        x_index = 1:2:length(R)*2;
        y_index = 2:2:length(R)*2;
        F_e(x_index) = F_e(x_index) +R*S(1,1)*Jmod;
        F_e(y_index) = F_e(y_index) +R*S(2,2)*Jmod;
    end
    idx = BoundaryBasis(eConn(e,:))';
    idx = global_id(:,idx);
    idx = idx(:); 
    F(idx) = F(idx) + F_e;
end
%% Dirichlet Boundary Condition
% There is one control point that is repeated. This means that we have to
% make their DOFs have the same displacements.
% This is done via penalty method
penalty = 1000;
P = Omega.get_point_cell;
PP = cell2mat(P(:));
repeated_points = intersect(find(PP(:,1) == -4),find(PP(:,2) == 4));
penalty_stiffness = penalty*[1 -1; -1 1];
penalty_ids = global_id(:,repeated_points);
K(penalty_ids(1,:),penalty_ids(1,:)) = ...
    K(penalty_ids(1,:),penalty_ids(1,:)) + penalty_stiffness;
K(penalty_ids(2,:),penalty_ids(2,:)) = ...
    K(penalty_ids(2,:),penalty_ids(2,:)) + penalty_stiffness;


% There are symmetry conditions to be enforced.
% The symmetry conditions say that:
% 1. There's no "flux" (the dot product of derivative with normal is zero)
% 2. The displacement normal to the plane of symmetry is zero
% Since zero derivativeis already natural in the formulation of the problem,
% We have to apply a Dirichlet BC on the symmetry side
% Let's identify the symmetric nodes
% There are symmetries in y = 0 and x = 4

Points = cell2mat(Points);
SymmetryBoundaries1 = find(abs(Points(:,1) - 0) < sqrt(eps)); % Careful with this logical
SymmetryBoundaries2 = find(abs(Points(:,2) - 0) < sqrt(eps)); % use abs(P - value) < eps
[global_id, ~] = BuildGlobalLocalMatrices(element_local_mapping,2);
SymmetryBoundaries1 = global_id(1,SymmetryBoundaries1); % x is symmetric
SymmetryBoundaries2 = global_id(2,SymmetryBoundaries2); % y is symmetric
SymmetryBoundaries = union(SymmetryBoundaries1,SymmetryBoundaries2);
SymmetryBoundaries = SymmetryBoundaries(:);
boundary_values = zeros(size(SymmetryBoundaries));
[d,F] = ApplyDirichletBCs(K,F,SymmetryBoundaries,boundary_values);

%% Post-Processing
dd = reshape(d,size(global_id));
scaling_factor = 100;
for i=1:length(Points)
Points(i,1:2) = Points(i,1:2) + scaling_factor*dd(:,i)';
end
PP = num2cell(Points,2);
PP = reshape(PP,size(Omega.PX));
DeformedOmega = Geometry('surf',Omega.pu,Omega.U,Omega.pv,Omega.V,PP);

% Stress Plotting

[~, NUMBER_OF_ELEMENTS] = size(element_local_mapping);
p = Omega.PolynomialOrder;
[quad_point_index, weights] = GenerateQuadPoints(p);
N_QUAD_POINTS = length(weights);
D = (YOUNG/(1-POISSON^2))*[1        POISSON     0; 
                           POISSON      1       0; 
                           0            0   (1-POISSON)/2];
Sigma = zeros(NUMBER_OF_ELEMENTS,N_QUAD_POINTS,3);
Coordinates = Sigma;
for e=1:NUMBER_OF_ELEMENTS
    u_range = element_ranges(e,:,1);
    v_range = element_ranges(e,:,2);
   for n=1:N_QUAD_POINTS
        IntegrationPoint = quad_point_index(n,:);
        [R, dR, J] = FastShape(Omega,IntegrationPoint,global_basis_index, ...
            element_local_mapping, element_ranges, e);
        qu = IntegrationPoint(1);
        qv = IntegrationPoint(2);
        u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2;        
        v = ((v_range(2)-v_range(1))*qv +(sum(v_range)))/2;    
        
        B = zeros(3,length(dR)*2);
        B(1,1:2:end) = dR(:,1);
        B(2,2:2:end) = dR(:,2);
        B(3,1:2:end) = dR(:,2);
        B(3,2:2:end) = dR(:,1);
        ActivePoints = element_local_mapping(:,e);
        ActiveBasis = global_id(:,ActivePoints);
        ActiveBasis = ActiveBasis(:);
        strain = B*d(ActiveBasis);
        stress = D*strain;
        Sigma(e,n,:) = stress(:);
        physical_coord = Omega.eval_point(u,v);
        Coordinates(e,n,1) = physical_coord.x;
        Coordinates(e,n,2) = physical_coord.y;
        Coordinates(e,n,3) = physical_coord.z;
    end
end
xx = Coordinates(:,:,1);
yy = Coordinates(:,:,2);
ss = Sigma(:,:,1);

figure(1)
scatter(xx(:),yy(:),10,ss(:))
% set(gca,'Color','k)
grid on
colormap(jet)
colorbar
        
        
        
        
        
        
        
        
        
        
        