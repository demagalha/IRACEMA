function arc = geo_arc(center,radius,plane,quadrant)

diameter = radius*2;

P = cell(1,3);
U = [0,0,0,1,1,1];

switch plane
    case 'xy'
        switch quadrant
            case 1        
                P{1} = [center(1) + diameter/2, center(2), center(3), 1];
                P{2} = [center(1) + diameter/2, center(2) + diameter/2, center(3), sqrt(2)/2];
                P{3} = [center(1), center(2) + diameter/2, center(3), 1];
            case 2
                P{1} = [center(1), center(2) + diameter/2, center(3), 1];
                P{2} = [center(1) - diameter/2, center(2) + diameter/2, center(3), sqrt(2)/2];
                P{3} = [center(1) - diameter/2, center(2), center(3), 1];
            case 3
                P{1} = [center(1) - diameter/2, center(2), center(3), 1];
                P{2} = [center(1) - diameter/2, center(2) - diameter/2, center(3), sqrt(2)/2];
                P{3} = [center(1), center(2) - diameter/2, center(3), 1];
            case 4
                P{1} = [center(1), center(2) - diameter/2, center(3), 1];
                P{2} = [center(1) + diameter/2, center(2) - diameter/2, center(3), sqrt(2)/2];
                P{3} = [center(1) + diameter/2, center(2), center(3), 1];
        end
end

arc = Geometry('curve',2,U,P);

end

