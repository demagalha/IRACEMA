function Model = BuildQuarterHemisphere(MidRadius,t)
a = 1/sqrt(2);
r = MidRadius-t/2;
R = MidRadius+t/2;
P{1} = [r 0 0 1]; P{2} = [R 0 0 1];
P{3} = [r 0 r a]; P{4} = [R 0 R a];
P{5} = [0 0 r 1]; P{6} = [0 0 R 1];
P{7} = [r r 0 a]; P{8} = [R R 0 a];
P{9} = [r r r 0.5]; P{10} = [R R R 0.5];
P{11} = [0 0 r a]; P{12} = [0 0 R a];
P{13} = [0 r 0 1]; P{14} = [0 R 0 1];
P{15} = [0 r r a]; P{16} = [0 R R a];
P{17} = [0 0 r 1]; P{18} = [0 0 R 1];

tmp = 1;
for i=1:3
    for j=1:3
        for k=1:2
            B{i,j,k} = P{tmp};
            tmp = tmp+1;
        end
    end
end
            
U = [0 0 0 1 1 1];
V = [0 0 0 1 1 1];
W = [0 0 1 1];
pu = 2;
pv = 2;
pw = 1;

Model = Geometry('volume',pu,U,pv,V,pw,W,B);
end