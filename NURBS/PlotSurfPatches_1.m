function [] = PlotSurfPatches_1(PX,PY,PZ,PW,nu,pu,U,nv,pv,V,render,a)

tam(1) = nu+1; tam(2) = nv+1;

Pw(1:tam(1),1:tam(2)) = CPOINT(0,0,0,0,1);

for i=1:tam(1)
    for j=1:tam(2)
        Pw(i,j) = CPOINT(PX(i,j),PY(i,j),PZ(i,j),PW(i,j),0);
    end
end

if strcmp(render,'coarse') == 1 
vv = linspace(V(1),V(end),50);

elseif strcmp(render,'fine') == 1
vv = linspace(V(1),V(end),200);

elseif strcmp(render,'medium') == 1
vv = linspace(V(1),V(end),100);
end

vv = [vv unique(V)];
vv = unique(sort(vv));

for j = numel(vv):-1:1
          S(j).x = 0;
          S(j).y = 0;
          S(j).z = 0;
end

for j=1:numel(vv)
    S(j) = SurfacePointRAT3(nu,pu,U,nv,pv,V,Pw,a,vv(j));
end

  
SX = zeros(1,numel(vv));
SY = zeros(1,numel(vv));
SZ = zeros(1,numel(vv));

for j=1:numel(vv)
    SX(j) = S(j).x;
    SY(j) = S(j).y;
    SZ(j) = S(j).z;
end
plot3(SX,SY,SZ,'-','color','black','LineWidth',2);

end