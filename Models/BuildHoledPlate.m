function Model = BuildHoledPlate(L,R)
P{1} = [0 0 0 1];
P{2} = [0 L 0 1];
P{3} = [L L 0 1];
C1 = Geometry('curve',1,[0 0 0.5 1 1],P);

B{1} = [L-R 0 0 1];
B{2} = [L-R R 0 sqrt(2)/2];
B{3} = [L R 0 1];
C2 = Geometry('curve',2,[0 0 0 1 1 1],B);

Model = geo_ruled(C1,C2,C1,C1);
end