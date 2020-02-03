function circle = geo_circle(center,diameter,plane)
 

P = cell(1,9);
U = [0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1];

switch plane
	case 'xy'
		P{1} = [center(1) + diameter/2, center(2), center(3), 1];
		P{2} = [center(1) + diameter/2, center(2) + diameter/2, center(3), sqrt(2)/2];
		P{3} = [center(1), center(2) + diameter/2, center(3), 1];
		P{4} = [center(1) - diameter/2, center(2) + diameter/2, center(3), sqrt(2)/2];
		P{5} = [center(1) - diameter/2, center(2), center(3), 1];
		P{6} = [center(1) - diameter/2, center(2) - diameter/2, center(3), sqrt(2)/2];
		P{7} = [center(1), center(2) - diameter/2, center(3), 1];
		P{8} = [center(1) + diameter/2, center(2) - diameter/2, center(3), sqrt(2)/2];
		P{9} = [center(1) + diameter/2, center(2), center(3), 1];
	case 'xz'
		P{1} = [center(1) + diameter/2, center(2), center(3), 1];
		P{2} = [center(1) + diameter/2, center(2), center(3) + diameter/2, sqrt(2)/2];
		P{3} = [center(1), center(2), center(3) + diameter/2, 1];
		P{4} = [center(1) - diameter/2, center(2)  center(3) + diameter/2, sqrt(2)/2];
		P{5} = [center(1) - diameter/2, center(2), center(3), 1];
		P{6} = [center(1) - diameter/2, center(2), center(3) - diameter/2, sqrt(2)/2];
		P{7} = [center(1), center(2), center(3) - diameter/2, 1];
		P{8} = [center(1) + diameter/2, center(2), center(3) - diameter/2, sqrt(2)/2];
		P{9} = [center(1) + diameter/2, center(2), center(3), 1];
	case 'yz'
		P{1} = [center(1), center(2) + diameter/2, center(3), 1];
		P{2} = [center(1), center(2) + diameter/2, center(3) + diameter/2, sqrt(2)/2];
		P{3} = [center(1), center(2), center(3) + diameter/2, 1];
		P{4} = [center(1), center(2) - diameter/2, center(3) + diameter/2, sqrt(2)/2];
		P{5} = [center(1), center(2) - diameter/2, center(3), 1];
		P{6} = [center(1), center(2) - diameter/2, center(3) - diameter/2, sqrt(2)/2];
		P{7} = [center(1), center(2), center(3) - diameter/2, 1];
		P{8} = [center(1), center(2) + diameter/2, center(3) - diameter/2, sqrt(2)/2];
		P{9} = [center(1), center(2) + diameter/2, center(3), 1];
	otherwise
		warning('invalid plane');
end

circle = Geometry('curve',2,U,P);

end
	