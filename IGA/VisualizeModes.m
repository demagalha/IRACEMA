function Modes = VisualizeModes(Model,autovector)
    [global_basis_index, element_local_mapping, element_ranges] = ...
        GetConnectivityArrays(Model);
    SOLUTION_DIMENSIONS = length(autovector)/length(global_basis_index);
   [ID, ~] = BuildGlobalLocalMatrices(element_local_mapping, ...
        SOLUTION_DIMENSIONS);   
    B = Model.get_point_cell;
    u = cell(size(B));
    comb = u;
    [sz1 sz2] = size(autovector);
    Modes = cell(sz2,1);
    pad = zeros(1,3-size(ID,1));
    for mode_n=1:sz2
        for i=1:size(ID,2)
            u{i} = [pad, autovector(ID(:,i),mode_n)', 0];
            comb{i} = B{i} +u{i};
        end
        switch Model.type
            case 'curve'
                Modes{mode_n} = Geometry('curve',Model.pu,Model.U,comb);
            case 'surf'
                Modes{mode_n} = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,comb);
            case 'volume'
                Modes{mode_n} = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
        end
    end    
end