function Model = BuildQuarterCilinder(OuterRadius,InnerRadius,Height)
R = OuterRadius;
r = InnerRadius;
h = Height;
a = 1/sqrt(2);

P{1} = [r 0 0 1]; P{2} = [R 0 0 1];
P{3} = [r 0 h 1]; P{4} = [R 0 h 1];
P{5} = [r r 0 a]; P{6} = [R R 0 a];
P{7} = [r r h a]; P{8} = [R R h a];
P{9} = [0 r 0 1]; P{10} = [0 R 0 1];
P{11} = [0 r h 1]; P{12} = [0 R h 1];
count = 1;

for i=1:3
    for k=1:2
        for j=1:2
            B{i,j,k} = P{count};
            count = count+1;
        end
    end
end
U = [0 0 0 1 1 1];
V = [0 0 1 1];
W = [0 0 1 1];

Model = Geometry('volume',2,U,1,V,1,W,B);
end