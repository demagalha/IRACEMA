clearvars
run startup.m
clc

%% Geometry Input
    a = 1;
    b = 1.5;
    % Ruled Surface Definition
        P1 = [0 0 0 1];
        P2 = [0 a 0 1];
        line1 = geo_line(P1,P2); % First line
        P3 = [b 0 0 1];
        P4 = [b a 0 1];
        line2 = geo_line(P3,P4); % Second line
        Omega = geo_ruled(line1,line2); % The Ruled Surface defined by lines 1 and 2
    % Subdomains Definition
        overlap = 0.75; % Percentage of intersection area
        P5 = [b*(1+overlap)/2 0 0 1];
        P6 = [b*(1+overlap)/2 a 0 1];
        line3 = geo_line(P5,P6);
        P7 = [b*(1-overlap)/2 0 0 1];
        P8 = [b*(1-overlap)/2 a 0 1];
        line4 = geo_line(P7,P8);
        Omega1 = geo_ruled(line1,line3); % Subdomain 1
        Omega2 = geo_ruled(line4,line2); % Subdomain 2
%% Refinement
    p = 2;
    Omega.DegreeElevate(p,1);
    Omega.DegreeElevate(p,2);
    Omega1.DegreeElevate(p,1);
    Omega1.DegreeElevate(p,2);
    Omega2.DegreeElevate(p,1);
    Omega2.DegreeElevate(p,2);
    h = 10;
    interval = linspace(0,1,h+2);
    interval = setdiff(interval,[0,1]);
%     Omega.KnotRefine(interval,1);
%     Omega.KnotRefine([1/3 2/3],2);
    Omega1.KnotRefine([.2 .4 .6 .8],1);
    Omega1.KnotRefine(interval,2);
    Omega2.KnotRefine([.2 .4 .6 .8],1);
    Omega2.KnotRefine(interval,2);
% We need to build inverse functions to communicate between the
% subdomains. We will also refine each domain right on the horizon between
% the overlapping part and the non overlapping part.

OmegaOneInverse = @(x,y) [y, 2*x/(b*(1+overlap))];
tmp = b*(1-overlap)/2;
OmegaTwoInverse = @(x,y) [y, (x-tmp)/(b-tmp)];

xx = Omega1.eval_point(1,1);
h2 = OmegaTwoInverse(xx.x,xx.y);
xx = Omega2.eval_point(0,0);
h1 = OmegaOneInverse(xx.x,xx.y);
Omega1.KnotRefine(h1(2),2);
Omega2.KnotRefine(h2(2),2);

%% Assembly
% Stiffness and Mass matrices
    NUMBER_OF_SOLUTION_DIMENTIONS = 1;
    [K, M] = PoissonProblemAssembly(Omega,NUMBER_OF_SOLUTION_DIMENTIONS);
    [K1,M1] = PoissonProblemAssembly(Omega1,NUMBER_OF_SOLUTION_DIMENTIONS);
    [K2,M2] = PoissonProblemAssembly(Omega2,NUMBER_OF_SOLUTION_DIMENTIONS);
% Boundary Definitions
    DomainBoundaries = GetBoundaries(Omega);
    SubdomainOneBoundaries = GetBoundaries(Omega1);
    SubdomainTwoBoundaries = GetBoundaries(Omega2);
     % Gamma 1 and 2 -> Boundaries defined by u = 0 and u = 1
        Gamma1 = DomainBoundaries{1};
        Gamma11 = SubdomainOneBoundaries{1};
        Gamma12 = SubdomainTwoBoundaries{1};
      
        Gamma2 = DomainBoundaries{2};
        Gamma21 = SubdomainOneBoundaries{2};
        Gamma22 = SubdomainTwoBoundaries{2};
    % Gammas 3 and 4 -> Boundaries defined by v = 0 and v = 1
        Gamma3 = DomainBoundaries{3};
        Gamma31 = SubdomainOneBoundaries{3};
        Gamma32 = SubdomainTwoBoundaries{3};
      
        Gamma4 = DomainBoundaries{4};
        Gamma41 = SubdomainOneBoundaries{4};
        Gamma42 = SubdomainTwoBoundaries{4};
    
%OBS: The overlapping boundaries in this case are:
%Omega1 -> Gamma41
%Omega2 -> Gamma32

% Boundary Restraints
    % Constrained Boundaries in this case are each boundary that is not
    % overlapping another domain
    % We have to run a unique() command so we don't constrain the same boundary
    % twice
    ConstrainedBoundaries = unique([Gamma1;Gamma2;Gamma3;Gamma4]);
    ConstrainedOneBoundaries = unique([Gamma11; Gamma21; Gamma31]);
    ConstrainedTwoBoundaries = unique([Gamma12; Gamma22; Gamma42]);

% Constraining the DOFs
    % First, let's "save" the original matrices
    KK = K; MM = M; 
    % Now we delete the constrained rows and columns from the matrices
    K(ConstrainedBoundaries,:) = [];
    K(:,ConstrainedBoundaries) = [];
    M(ConstrainedBoundaries,:) = [];
    M(:,ConstrainedBoundaries) = [];
    
    % Repeating for the SubDomains
    KK1 = K1; MM1 = M1;
    K1(ConstrainedOneBoundaries,:) = [];
    K1(:,ConstrainedOneBoundaries) = [];
    M1(ConstrainedOneBoundaries,:) = [];
    M1(:,ConstrainedOneBoundaries) = [];

    KK2 = K2; MM2 = M2;
    K2(ConstrainedTwoBoundaries,:) = [];
    K2(:,ConstrainedTwoBoundaries) = [];
    M2(ConstrainedTwoBoundaries,:) = [];
    M2(:,ConstrainedTwoBoundaries) = [];

% Solving the Eigenvalue Problems
    K = sparse(K);
    M = sparse(M);
    [EigenVector, EigenValues] = eigs(K,M,1,'sm');
    
    K1 = sparse(K1);
    M1 = sparse(M1);
    [EigenVector1, EigenValues1] = eigs(K1,M1,1,'sm');
    K1 = KK1;
    
    K2 = sparse(K2);
    M2 = sparse(M2);
    [EigenVector2, EigenValues2] = eigs(K2,M2,1,'sm');
    K2 = KK2;
%% Post-Processing
% We have to re-add the constrained DOFs in the right position as 0s in the
% EigenVector. Fortunately, we have made a function for that.
EigenVector = BoundariesPostProcess(EigenVector,ConstrainedBoundaries);
EigenVector1 = BoundariesPostProcess(EigenVector1,ConstrainedOneBoundaries);
EigenVector2 = BoundariesPostProcess(EigenVector2,ConstrainedTwoBoundaries);

% Visualization
    OmegaVisualizationCell = VisualizeModes(Omega,EigenVector);
    Omega1VisualizationCell = VisualizeModes(Omega1,EigenVector1);
    Omega2VisualizationCell = VisualizeModes(Omega2,EigenVector2);

%% Schwarz Algorithm Loops
% With the inital problem setted up, we now cycle through a series of Lui
% Boundary Conditions in Omega1 and Omega2 until we satisfy a convergence
% condition
L2norm = 200;
count = 0;
while (L2norm < 1e-3 || count < 15)
    % Lui Boundary Condition for Omega 1
    K1 = LuiBC_2D(Omega1VisualizationCell{1},Omega2VisualizationCell{1}, ...
        KK1,Gamma41,4,OmegaTwoInverse);
        M1 = MM1;
        K1(ConstrainedOneBoundaries,:) = [];
        K1(:,ConstrainedOneBoundaries) = [];
        M1(ConstrainedOneBoundaries,:) = [];
        M1(:,ConstrainedOneBoundaries) = [];
        K1 = sparse(K1);
        M1 = sparse(M1);
        [EigenVector1, EigenValues1] = eigs(K1,M1,1,'sm');
        if abs(min(EigenVector1)) > abs(max(EigenVector1))
            EigenVector1 = -EigenVector1;
        end
        K1 = KK1;
    EigenVector1 = BoundariesPostProcess(EigenVector1,ConstrainedOneBoundaries);
    Omega1VisualizationCell = VisualizeModes(Omega1,EigenVector1);
    
    % Lui Boundary Condition for Omega 2
    K2 = LuiBC_2D(Omega2VisualizationCell{1},Omega1VisualizationCell{1}, ...
        KK2,Gamma32,3,OmegaOneInverse);
        M2 = MM2;
        K2(ConstrainedTwoBoundaries,:) = [];
        K2(:,ConstrainedTwoBoundaries) = [];
        M2(ConstrainedTwoBoundaries,:) = [];
        M2(:,ConstrainedTwoBoundaries) = [];
        K2 = sparse(K2);
        M2 = sparse(M2);
        [EigenVector2, EigenValues2] = eigs(K2,M2,1,'sm');
        if abs(min(EigenVector2)) > abs(max(EigenVector2))
            EigenVector2 = -EigenVector2;
        end
        K2 = KK2;
    EigenVector2 = BoundariesPostProcess(EigenVector2,ConstrainedTwoBoundaries);
    Omega2VisualizationCell = VisualizeModes(Omega2,EigenVector2);
 
    % Calculation of the L2 norm of the difference
    % We have to find which elements of Omega1 lie in the overlapping
    % domain (elements of Omega2 aren't necessary)
    % We define these elements as the elements which the support start at
    % the overlap
    % So, to find which elements belong to the overlapping domain, we have
    % first to find the basis functions that end/start at the domain
    % separation
    b1 = Omega2.eval_point(0,0);
    overlap_start = OmegaOneInverse(b1.x,b1.y);
    % We know the overlap ends at u = 1, v = 1 in Omega1
    U = Omega1.U;
    V = Omega1.V;
    pu = Omega1.pu;
    pv = Omega1.pv;
    uspan = FindSpanLinear(length(U)-pu-1,pu,overlap_start(1),U);
    u_start = uspan+1-pu;
    vspan = FindSpanLinear(length(V)-pv-1,pv,overlap_start(2),V);
    v_start = vspan+1-pv;
    [global_basis_index, element_local_mapping, element_ranges] = ...
         GetConnectivityArrays(Omega1);
     bool1 = global_basis_index(:,1) >= u_start;
     bool2 = global_basis_index(:,2) >= v_start;
     bool = bool1 & bool2;
     basis_start = find(bool==1);
     [~, elements, ~] = intersect(element_local_mapping(1,:),basis_start);
     
     p = Omega1.PolynomialOrder;
     L2norm = 0;
     [quad_point_index, weights] = GenerateQuadPoints(p);
     for ee=1:length(elements)
         e = elements(ee);
         for n=1:length(weights)
             qu = quad_point_index(n,1);
             qv = quad_point_index(n,2);
             support = global_basis_index(element_local_mapping(:,e),:);
             u_range = element_ranges(e,:,1);
             v_range = element_ranges(e,:,2);
             u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric
             v = ((v_range(2)-v_range(1))*qv +(sum(v_range)))/2;
             zz = Omega1VisualizationCell{1}.eval_point(u,v);
             z1 = zz.z;
             parametric_2 = OmegaTwoInverse(zz.x,zz.y);
             zz = Omega2VisualizationCell{1}.eval_point(parametric_2(1),parametric_2(2));
             z2 = zz.z;
             L2norm = L2norm + weights(n)*(z1 - z2)^2;
         end
     end
     L2norm = sqrt(L2norm)
     count = count+1
end
figure(1)
Omega1VisualizationCell{1}.plot_geo;
Omega2VisualizationCell{1}.plot_geo;
alpha(0.6)