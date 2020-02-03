function Model = BuildCilinder(OuterRadius,InnerRadius,Height)
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
P{13} = [-r r 0 a]; P{14} = [-R R 0 a];
P{15} = [-r r h a]; P{16} = [-R R h a];
P{17} = [-r 0 0 1]; P{18} = [-R 0 0 1];
P{19} = [-r 0 h 1]; P{20} = [-R 0 h 1];
P{21} = [-r -r 0 a]; P{22} = [-R -R 0 a];
P{23} = [-r -r h a]; P{24} = [-R -R h a];
P{25} = [0 -r 0 1]; P{26} = [0 -R 0 1];
P{27} = [0 -r h 1]; P{28} = [0 -R h 1];
P{29} = [r -r 0 a]; P{30} = [R -R 0 a];
P{31} = [r -r h a]; P{32} = [R -R h a];
P{33} = [r 0 0 1]; P{34} = [R 0 0 1];
P{35} = [r 0 h 1]; P{36} = [R 0 h 1];

count = 1;

for i=1:9
    for k=1:2
        for j=1:2
            B{i,j,k} = P{count};
            count = count+1;
        end
    end
end
U = [0 0 0 1 1 2 2 3 3 4 4 4];
V = [0 0 1 1];
W = [0 0 1 1];

Model = Geometry('volume',2,U,1,V,1,W,B);
end