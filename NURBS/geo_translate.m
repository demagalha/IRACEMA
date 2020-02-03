function translate = geo_translate(Model,vector)

P = Model.get_point_cell;
tam = size(P);
vec = [vector 0];
switch Model.type
    case 'volume'
        for i=1:tam(1)
            for j=1:tam(2)
                for k=1:tam(3)
                    P{i,j,k} = P{i,j,k} + vec;
                end
            end
        end
        translate = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,P);
end

end