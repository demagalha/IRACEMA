function scaling = geo_scaling(Model,vector)

S = [vector(1) 0 0; 0 vector(2) 0; 0 0 vector(3)];

P = Model.get_point_cell;
tam = size(P);

switch Model.type
    case 'surf'
        for i=1:tam(1)
            for j=1:tam(2)
                B{i,j} = [P{i,j}(1); P{i,j}(2); P{i,j}(3)];
                B{i,j} = S*B{i,j};
                B{i,j}(4) = P{i,j}(4);
            end
        end
        scaling = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,B);
        
    case 'curve'
        for i=1:numel(P)
            B{i} = [P{i}(1);P{i}(2);P{i}(3)];
            B{i} = S*B{i};
            B{i}(4) = P{i}(4);
        end
        scaling = Geometry('curve',Model.pu,Model.U,B);
end

end
