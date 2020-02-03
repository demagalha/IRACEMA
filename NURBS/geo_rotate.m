function rotate = geo_rotate(Model,eixo,theta)

P = Model.get_point_cell;
tam = size(P);
rot = [cos(theta) + eixo(1)*eixo(1)*(1-cos(theta)), eixo(1)*eixo(2)*(1-cos(theta)) - eixo(3)*sin(theta), eixo(1)*eixo(3)*(1-cos(theta)) + eixo(2)*sin(theta);
                   eixo(2)*eixo(1)*(1-cos(theta)) + eixo(3)*sin(theta), cos(theta) + eixo(2)*eixo(2)*(1-cos(theta)), eixo(2)*eixo(3)*(1-cos(theta)) - eixo(1)*sin(theta);
                   eixo(3)*eixo(1)*(1-cos(theta)) - eixo(2)*sin(theta), eixo(3)*eixo(2)*(1-cos(theta)) + eixo(1)*sin(theta), cos(theta) + eixo(3)*eixo(3)*(1-cos(theta))];
switch Model.type
    case 'surf'
        for i=1:tam(1)
            for j=1:tam(2)
                B{i,j} = [P{i,j}(1); P{i,j}(2); P{i,j}(3)];
                B{i,j} = rot*B{i,j};
                B{i,j}(4) = P{i,j}(4);
            end
        end
        rotate = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,B);
        
    case 'curve'
        for i=1:numel(P)
            B{i} = [P{i}(1);P{i}(2);P{i}(3)];
            B{i} = rot*B{i};
            B{i}(4) = P{i}(4);
        end
        rotate = Geometry('curve',Model.pu,Model.U,B);
        
    case 'volume'
        for i=1:tam(1)
            for j=1:tam(2)
                for k=1:tam(3)
                    B{i,j,k} = [P{i,j,k}(1);P{i,j,k}(2);P{i,j,k}(3)];
                    B{i,j,k} = rot*B{i,j,k};
                    B{i,j,k}(4) = P{i,j,k}(4);
                end
            end
        end
        rotate = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,B);
end

end
        
        