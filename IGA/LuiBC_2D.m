function K = LuiBC_2D(GeometryObj1,GeometryObj2,K,BoundaryBasis,Boundary, ...
    InverseFunction)
    [global_basis_index, element_local_mapping, element_ranges] = ...
         GetConnectivityArrays(GeometryObj1);
    directions = [1 2];
    QUAD_DIRECTION = setdiff(directions,round(Boundary/2));
    BOUNDARY_DIRECTION = round(Boundary/2);
    
    p = GeometryObj1.PolynomialOrder;
    Knots = GeometryObj1.KnotVectorCell;
    
    pq = p(QUAD_DIRECTION);
    Knot = Knots{QUAD_DIRECTION};
    [boundary_ranges, eConn] = KnotConnectivity(pq,Knot);
    [NUMBER_OF_ELEMENTS, ELEMENT_DOFS] = size(eConn);
    
    [q, w] = getGP(pq);
    Points = GeometryObj1.get_point_cell;
    BPoints = Points(:);
    BPoints = BPoints(BoundaryBasis); 
    BPoints = BPoints';
    BoundaryObj = Geometry('curve',pq,Knot,BPoints);
%     BPoints = cell2mat(BPoints);
    for e=1:NUMBER_OF_ELEMENTS
        K_e = zeros(ELEMENT_DOFS,1);
        B_RANGE = boundary_ranges(e,:);
        parametric_coordinates = zeros(2,1);
        for i=1:length(q)
            IntegrationPoint = q(i);
            [R, dR, J] = FastBoundaryShape(GeometryObj1,IntegrationPoint, ...
                QUAD_DIRECTION, B_RANGE, Points, e, eConn);
            u1 = ((B_RANGE(2)-B_RANGE(1))*q(i) +(sum(B_RANGE)))/2;
            v1 = 1-mod(Boundary,2);
            parametric_coordinates(QUAD_DIRECTION) = u1;
            parametric_coordinates(BOUNDARY_DIRECTION) = v1;
            u1 = parametric_coordinates(1);
            v1 = parametric_coordinates(2);
            normal = CurveNurbsNormal(BoundaryObj, ...
                parametric_coordinates(QUAD_DIRECTION));
            SpatialCoordinates = GeometryObj1.eval_point(u1,v1);
            x = SpatialCoordinates.x;
            y = SpatialCoordinates.y;
            parametric = InverseFunction(x,y);
            u2 = parametric(1);
            v2 = parametric(2);
            SubDomainSpatial = GeometryObj2.eval_point(u2,v2);
            z2 = SubDomainSpatial.z;
            g_ell = z2;
            
            hx = sqrt(eps)*x;
            hy = sqrt(eps)*y;
            perturbation = InverseFunction(x+hx,y+hy);
            SubDomainPerturbed = GeometryObj2.eval_point(perturbation(1), ...
                perturbation(2));
            dx = SubDomainPerturbed.x - x;
            dy = SubDomainPerturbed.y - y;
            dz = SubDomainPerturbed.z - z2;
            dzdx = dz/dx;
            dzdy = dz/dy;
            du2 = [dzdx; dzdy; 0];
            f_ell = dot(du2,normal);
            K_e = K_e +w(i)*J*(f_ell*(R*R') -(1-g_ell)*(dR*R'));
        end
        idx = BoundaryBasis(eConn(e,:));
        K(idx,idx) = K(idx,idx) + K_e;
    end         
end