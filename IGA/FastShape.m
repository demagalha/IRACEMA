function [R, dR, J] = FastShape(GeometryObject,IntegrationPoint, ...
    global_basis_index, element_local_mapping,element)
    switch GeometryObject.type
        case 'curve'
            [R, dR, J] = FastShape1D(GeometryObject,IntegrationPoint, ... 
                global_basis_index, element_local_mapping,element);
        case 'surface'
            [R, dR, J] = FastShape2D(GeometryObject,IntegrationPoint, ... 
                global_basis_index, element_local_mapping,element);
        case 'volume'
            [R, dR, J] = FastShape3D(GeometryOject,IntegrationPoint, ...
                global_basis_index, element_local_mapping,element);
    end
end