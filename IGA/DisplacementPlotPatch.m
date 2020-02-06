function [DX, DY, DZ] = DisplacementPlotPatch(autovector,ID,Model,option,render)

ModelDisplacement = DisplacementModel(autovector,ID,Model); 

%extract faces from volumetric model (geometry) and (displacement as geometry)
SurfacesModel = ExtractFaces(Model);
SurfacesDisplacement = ExtractFaces(ModelDisplacement);

S = cell(6,1);
D = cell(6,1);

SX = cell(6,1);
SY = cell(6,1);
SZ = cell(6,1);
DX = cell(6,1);
DY = cell(6,1);
DZ = cell(6,1);
DT = cell(6,1);


S{1} = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,SurfacesModel{1}.P);
S{6} = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,SurfacesModel{6}.P);
D{1} = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,SurfacesDisplacement{1}.P);
D{6} = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,SurfacesDisplacement{6}.P);

S{4} = Geometry('surf',Model.pw,Model.W,Model.pv,Model.V,SurfacesModel{4}.P);
S{3} = Geometry('surf',Model.pw,Model.W,Model.pv,Model.V,SurfacesModel{3}.P);
D{4} = Geometry('surf',Model.pw,Model.W,Model.pv,Model.V,SurfacesDisplacement{4}.P);
D{3} = Geometry('surf',Model.pw,Model.W,Model.pv,Model.V,SurfacesDisplacement{3}.P);

S{2} = Geometry('surf',Model.pu,Model.U,Model.pw,Model.W,SurfacesModel{2}.P);
S{5} = Geometry('surf',Model.pu,Model.U,Model.pw,Model.W,SurfacesModel{5}.P);
D{2} = Geometry('surf',Model.pu,Model.U,Model.pw,Model.W,SurfacesDisplacement{2}.P);
D{5} = Geometry('surf',Model.pu,Model.U,Model.pw,Model.W,SurfacesDisplacement{5}.P);


max_points = 0;
switch render
    case 'coarse'
        max_points = 50;
    case 'medium'
        max_points = 100;
    case 'fine'
        max_points = 200;
    otherwise
        disp('invalid input for discretization');
end


for k=1:6
    
    uu = linspace(D{k}.U(1),D{k}.U(end),max_points);
    vv = linspace(D{k}.V(1),D{k}.V(end),max_points);
    
    uu = [uu unique(D{k}.U)];
    uu = unique(sort(uu));
    vv = [vv unique(D{k}.V)];
    vv = unique(sort(vv));
    
    SX{k} = zeros(numel(uu),numel(vv));
    SY{k} = SX{k};
    SZ{k} = SX{k};
    
    DX{k} = SX{k};
    DY{k} = SX{k};
    DZ{k} = SX{k};
    DT{k} = SX{k};
      
    for i=1:numel(uu)
        for j=1:numel(vv)
        [SX{k}(i,j),SY{k}(i,j),SZ{k}(i,j)] = SurfacePointRAT(S{k},uu(i),vv(j));
        [DX{k}(i,j),DY{k}(i,j),DZ{k}(i,j)] = SurfacePointRAT(D{k},uu(i),vv(j));
        DT{k}(i,j) = sqrt((DX{k}(i,j))^2 + (DY{k}(i,j))^2 + (DZ{k}(i,j))^2);
        end
    end
    
end

switch option
    case 'dx'
        minvalue = min([DX{1}(:);DX{2}(:);DX{3}(:);DX{4}(:);DX{5}(:);DX{6}(:)]);
        maxvalue = max([DX{1}(:);DX{2}(:);DX{3}(:);DX{4}(:);DX{5}(:);DX{6}(:)]);
        
        surf(SX{1},SY{1},SZ{1},DX{1},'edgecolor','none','FaceColor','interp');
        hold all;
        for i=2:6
            surf(SX{i},SY{i},SZ{i},DX{i},'edgecolor','none','FaceColor','interp');
        end
        
    case 'dy'
        minvalue = min([DY{1}(:);DY{2}(:);DY{3}(:);DY{4}(:);DY{5}(:);DY{6}(:)]);
        maxvalue = max([DY{1}(:);DY{2}(:);DY{3}(:);DY{4}(:);DY{5}(:);DY{6}(:)]);
        
        surf(SX{1},SY{1},SZ{1},DY{1},'edgecolor','none','FaceColor','interp');
        hold all;
        for i=2:6
            surf(SX{i},SY{i},SZ{i},DY{i},'edgecolor','none','FaceColor','interp');
        end
        
    case 'dz'
        minvalue = min([DZ{1}(:);DZ{2}(:);DZ{3}(:);DZ{4}(:);DZ{5}(:);DZ{6}(:)]);
        maxvalue = max([DZ{1}(:);DZ{2}(:);DZ{3}(:);DZ{4}(:);DZ{5}(:);DZ{6}(:)]);
        
        surf(SX{1},SY{1},SZ{1},DZ{1},'edgecolor','none','FaceColor','interp');
        hold all;
        for i=2:6
            surf(SX{i},SY{i},SZ{i},DZ{i},'edgecolor','none','FaceColor','interp');
        end
        
    case 'd'
        minvalue = min([DT{1}(:);DT{2}(:);DT{3}(:);DT{4}(:);DT{5}(:);DT{6}(:)]);
        maxvalue = max([DT{1}(:);DT{2}(:);DT{3}(:);DT{4}(:);DT{5}(:);DT{6}(:)]);
        
        surf(SX{1},SY{1},SZ{1},DT{1},'edgecolor','none','FaceColor','interp');
        hold all;
        for i=2:6
            surf(SX{i},SY{i},SZ{i},DT{i},'edgecolor','none','FaceColor','interp');
        end
        
    otherwise
        disp('Invalid input');
end
 
caxis([minvalue,maxvalue]);
colorbar;
colormap(jet);


end