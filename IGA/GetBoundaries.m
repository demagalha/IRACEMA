function Boundaries = GetBoundaries(GeometryObj)
    [global_basis_index, element_local_mapping, ~] = ...
        GetConnectivityArrays(GeometryObj);
    [~, sz] = size(global_basis_index);
    Boundaries = cell(sz*2,2);
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
        Boundaries{g,1} = basis_funs;
        [~, idx, ~] = intersect(element_local_mapping,basis_funs);
        [~, boundary_elements] = ind2sub(size(element_local_mapping),idx);
        boundary_elements = unique(boundary_elements);
        Boundaries{g,2} = boundary_elements;
        g = g+1;
    end
end