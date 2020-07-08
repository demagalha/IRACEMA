function [global_basis_index, element_local_mapping, global_id, local_matrix] = GetConnectivityArrays(GeometryObj)
    switch GeometryObj.type
        case 'curve'
            [global_basis_index, element_local_mapping] = CurveConnectivity(GeometryObj);
        case 'surface'
            [global_basis_index, element_local_mapping]  = SurfaceConnectivity(GeometryObj);
        case 'volume'
            [global_basis_index, element_local_mapping]  = VolumeConnectivity(GeometryObj);
    end
end