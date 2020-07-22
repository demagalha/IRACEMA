function Boundaries = GetBoundaries(GeometryObj)
    [global_basis_index, element_local_mapping, ~] = ...
        GetConnectivityArrays(GeometryObj);
    [~, sz] = size(global_basis_index);
    Boundaries = cell(sz*2,1);
    g = 1;
    basis_max = max(global_basis_index);
    while g <= sz*2
        column = round(g/2);
        if mod(g,2)
            basis = 1;
        else
            basis = basis_max(column);
        end
        basis_funs = find(global_basis_index(:,column) == basis);
        Boundaries{g} = basis_funs;
        [~, idx, ~] = intersect(element_local_mapping,basis_funs);
        g = g+1;
    end
end