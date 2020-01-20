function [] = PlotSurfPatches_2(PX,PY,PZ,PW,nu,pu,U,nv,pv,V,render,a,PlotStruct)

tam(1) = nu+1; tam(2) = nv+1;

Pw(1:tam(1),1:tam(2)) = CPOINT(0,0,0,0,1);

for i=1:tam(1)
    for j=1:tam(2)
        Pw(i,j) = CPOINT(PX(i,j),PY(i,j),PZ(i,j),PW(i,j),0);
    end
end

if strcmp(render,'coarse') == 1
    
uu = linspace(U(1),U(end),50);

elseif strcmp(render,'fine') == 1
uu = linspace(U(1),U(end),200);

elseif strcmp(render,'medium') == 1
uu = linspace(U(1),U(end),100);
end

uu = [uu unique(U)];
uu = unique(sort(uu));

for i = numel(uu):-1:1
    S(i).x = 0;
    S(i).y = 0;
    S(i).z = 0;
end
  
 for i=1:numel(uu)
    S(i) = SurfacePointRAT3(nu,pu,U,nv,pv,V,Pw,uu(i),a);
 end
  
SX = zeros(numel(uu),1);
SY = zeros(numel(uu),1);
SZ = zeros(numel(uu),1);

for i=1:numel(uu)
    SX(i) = S(i).x;
    SY(i) = S(i).y;
    SZ(i) = S(i).z;
end


plot3(SX,SY,SZ,'-','color',PlotStruct.IsoRGB,'LineWidth',PlotStruct.LineSize);

end