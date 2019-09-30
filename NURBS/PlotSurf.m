function [] = PlotSurf(PX,PY,PZ,PW,nu,pu,U,nv,pv,V,render)

tam(1) = nu+1; tam(2) = nv+1;

Pw(1:tam(1),1:tam(2)) = CPOINT(0,0,0,0,1);

for i=1:tam(1)
    for j=1:tam(2)
        Pw(i,j) = CPOINT(PX(i,j),PY(i,j),PZ(i,j),PW(i,j),0);
    end
end

if strcmp(render,'coarse') == 1
    
uu = linspace(U(1),U(end),50);
vv = linspace(V(1),V(end),50);

elseif strcmp(render,'fine') == 1
uu = linspace(U(1),U(end),200);
vv = linspace(V(1),V(end),200);

elseif strcmp(render,'medium') == 1
uu = linspace(U(1),U(end),100);
vv = linspace(V(1),V(end),100);
end

uu = [uu unique(U)];
uu = unique(sort(uu));

vv = [vv unique(V)];
vv = unique(sort(vv));

 for i = numel(uu):-1:1
      for j = numel(vv):-1:1
          S(i,j).x = 0;
          S(i,j).y = 0;
          S(i,j).z = 0;
      end
 end
  
 for i=1:numel(uu)
      for j=1:numel(vv)
      S(i,j) = SurfacePointRAT3(nu,pu,U,nv,pv,V,Pw,uu(i),vv(j));
      end
 end
  
 SX = zeros(numel(uu),numel(vv));
 SY = zeros(numel(uu),numel(vv));
 SZ = zeros(numel(uu),numel(vv));
for i=1:numel(uu)
    for j=1:numel(vv)
        SX(i,j) = S(i,j).x;
        SY(i,j) = S(i,j).y;
        SZ(i,j) = S(i,j).z;
    end
end


h = surf(SX,SY,SZ);
set(h,'edgecolor','none','facecolor',[.1 .9 .1])


end