function line = geo_line(point1,point2)

U = [0 0 1 1];
P{1} = [point1 1];
P{2} = [point2 1];

line = Geometry('curve',1,U,P);

end