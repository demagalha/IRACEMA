function [sigx] = PlotStress2D(Model,ModelDisplacement,YOUNG,POISSON,type)


max_points = 10;
uu = linspace(Model.U(1),Model.U(end),max_points);
vv = linspace(Model.V(1),Model.V(end),max_points);
uu = [uu unique(Model.U)];
uu = unique(sort(uu));
vv = [vv unique(Model.V)];
vv = unique(sort(vv));

SX = zeros(numel(uu),numel(vv));
SY = zeros(numel(uu),numel(vv));
SZ = SY;
sigx = SY;
sigy = sigx;
sigxy = sigx;

DX = ModelDisplacement.PX;
DY = ModelDisplacement.PY;
P = Model.get_point_cell;

for i=1:numel(uu)
    for j=1:numel(vv)
        [SX(i,j),SY(i,j),SZ(i,j)] = SurfacePointRAT(Model,uu(i),vv(j));
        sigma = CalculateStressVector2d(P,Model,DX,DY,uu(i),vv(j),YOUNG,POISSON);
%         if isnan(sigma(1))
%             sigma = CalculateStressVector2d(P,Model,DX,DY,uu(i)-0.01,vv(j),YOUNG,POISSON);
%             uu(i)
%             vv(j)
%         end
        sigx(i,j) = abs(sigma(1));
        sigy(i,j) = abs(sigma(2));
        sigxy(i,j) = abs(sigma(3));
    end
end

switch type
    case 'sigx'     
        surf(SX,SY,SZ,sigx,'edgecolor','none','FaceColor','interp');
    case 'sigy'
        surf(SX,SY,SZ,sigy,'edgecolor','none','FaceColor','interp');
    case 'sigxy'
        surf(SX,SY,SZ,sigxy,'edgecolor','none','FaceColor','interp');
end

colorbar;
colormap(jet);
 
end