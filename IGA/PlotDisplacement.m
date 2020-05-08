function PlotDisplacementContour(autovector,ID,Model)

T = VisualizeModes(Model,autovector,ID);
T = T{1};
max_points = 100;
uu = linspace(Model.U(1),Model.U(end),max_points);
vv = linspace(Model.V(1),Model.V(end),max_points);
uu = [uu unique(Model.U)];
uu = unique(sort(uu));
vv = [vv unique(Model.V)];
vv = unique(sort(vv));

SX = zeros(numel(uu),numel(vv)); SY = SX; SZ = SX; T_ = SX;

for i=1:numel(uu)
	for j=1:numel(vv)
		[SX(i,j),SY(i,j),SZ(i,j)] = SurfacePointRAT(Model,uu(i),vv(j));
		[~, ~, T_(i,j)] = SurfacePointRAT(T,uu(i),vv(j));
	end
end

surf(SX,SY,SZ,T_,'edgecolor','none','FaceColor','interp');
colorbar;
colormap(jet);