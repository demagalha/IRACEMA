function B = geo_rotate_points(P,eixo,theta,type)


tam = size(P);
rot = [cos(theta) + eixo(1)*eixo(1)*(1-cos(theta)), eixo(1)*eixo(2)*(1-cos(theta)) - eixo(3)*sin(theta), eixo(1)*eixo(3)*(1-cos(theta)) + eixo(2)*sin(theta);
                   eixo(2)*eixo(1)*(1-cos(theta)) + eixo(3)*sin(theta), cos(theta) + eixo(2)*eixo(2)*(1-cos(theta)), eixo(2)*eixo(3)*(1-cos(theta)) - eixo(1)*sin(theta);
                   eixo(3)*eixo(1)*(1-cos(theta)) - eixo(2)*sin(theta), eixo(3)*eixo(2)*(1-cos(theta)) + eixo(1)*sin(theta), cos(theta) + eixo(3)*eixo(3)*(1-cos(theta))];
switch type
    case 'surf'
        for i=1:tam(1)
            for j=1:tam(2)
                B{i,j} = [P{i,j}(1); P{i,j}(2); P{i,j}(3)];
                B{i,j} = rot*B{i,j};
                B{i,j}(4) = P{i,j}(4);
            end
        end
        
    case 'curve'
        for i=1:numel(P)
            B{i} = [P{i}(1);P{i}(2);P{i}(3)];
            B{i} = rot*B{i};
            B{i}(4) = P{i}(4);
        end
end

end
        
        