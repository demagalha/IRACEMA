function F = ApplyNeumannBCs(GeometryObj,K,F,BoundaryBasis,Boundary, ...
    boundary_values)
    switch GeometryObj.type
        case 'curve'
            [K,F] = NeumannBCs_1D(GeometryObj,K,F,BoundaryBasis,Boundary, ...
                boundary_values);
        case 'surf'
            [K,F] = NeumannBCs_2D(GeometryObj,K,F,BoundaryBasis,Boundary, ...
                boundary_values);
        case 'volume'
            [K,F] = NeumannBCs_3D(GeometryObj,K,F,BoundaryBasis,Boundary, ...
                boundary_values);
    end
end