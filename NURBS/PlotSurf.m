function [] = PlotSurf(Model,render,Urange,Vrange)

switch nargin
    case 2
        Urange = [Model.U(1) Model.U(end)];
        Vrange = [Model.V(1) Model.V(end)];
    case 3
        Vrange = [Model.V(1) Model.V(end)];
        if length(Urange) ~= 2
            error('Urange format is not right, should be [Umin Umax]')
        end
    case 4
        if length(Urange) ~= 2
            error('Urange format is not right, should be [Umin Umax]');
        end
        if length(Vrange) ~= 2
            error('Vrange format is not right, should be [Vmin Vmax]');
        end
end

switch render
    case 'coarse'    
        uu = linspace(Urange(1),Urange(2),50);
        vv = linspace(Vrange(1),Vrange(2),50);

    case 'fine'
        uu = linspace(Urange(1),Urange(2),200);
        vv = linspace(Vrange(1),Vrange(2),200);
        
    case 'medium'
        uu = linspace(Urange(1),Urange(2),100);
        vv = linspace(Vrange(1),Vrange(2),100);
end


if Urange(1) == Model.U(1) && Urange(2) == Model.U(end)
    uu = [uu unique(Model.U)];
    uu = unique(sort(uu));
end

if Vrange(1) == Model.V(1) && Vrange(2) == Model.V(end)
    vv = [vv unique(Model.V)];
    vv = unique(sort(vv));
end

  
SX = zeros(numel(uu),numel(vv));
SY = zeros(numel(uu),numel(vv));
SZ = zeros(numel(uu),numel(vv));

for i=1:numel(uu)
    for j=1:numel(vv)
        [SX(i,j),SY(i,j),SZ(i,j)] = SurfacePointRAT(Model,uu(i),vv(j));
    end
end


h = surf(SX,SY,SZ);
set(h,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong')


end