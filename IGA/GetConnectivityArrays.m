function [global_basis_index, element_local_mapping, element_ranges] = ...
    GetConnectivityArrays(GeometryObj)
    switch GeometryObj.type
        case 'curve'
            [global_basis_index, element_local_mapping, element_ranges] ...
                = CurveConnectivity(GeometryObj);
        case 'surf'
            [global_basis_index, element_local_mapping, element_ranges] ...
                = SurfaceConnectivity(GeometryObj);
        case 'volume'
            [global_basis_index, element_local_mapping, element_ranges] ...
                = VolumeConnectivity(GeometryObj);
    end
end