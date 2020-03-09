function [] = PlotRatVolRange(Model,Urange,Vrange,Wrange,render,cpoints,isolines)

switch render
    case 'coarse'
        points = 50;
    case 'medium'
        points = 100;
    case 'fine'
        points = 200;
end


PX = Model.PX; PY = Model.PY; PZ = Model.PZ;

Pwx = Model.PX .* Model.weight;
Pwy = Model.PY .* Model.weight;
Pwz = Model.PZ .* Model.weight;

u = linspace(Urange(1),Urange(2),points);
v = linspace(Vrange(1),Vrange(2),points);
w = linspace(Wrange(1),Wrange(2),points);

 %%%%%%%%%%%%%%%%%%% face 1
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 
SX = zeros(points,points);
SY = zeros(points,points);
SZ = zeros(points,points);
 
for i=1:points
    for j=1:points
         [SX(i,j),SY(i,j),SZ(i,j)] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u(i),v(j),w(1));
    end
end
 
surf(SX,SY,SZ,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong');
 
hold on;
 
clear SX SY SZ
 
  %%%%%%%%%%%%%%%%%%% face 6
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v         w=end |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 
SX = zeros(points,points);
SY = zeros(points,points);
SZ = zeros(points,points);
 
for i=1:points
    for j=1:points
        [SX(i,j),SY(i,j),SZ(i,j)] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u(i),v(j),w(end));
    end
end
 
surf(SX,SY,SZ,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong');

clear SX SY SZ

 %%%%%%%%%%%%%%%%%%% face 4
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->w----
 
SX = zeros(points);
SY = zeros(points);
SZ = zeros(points);

for i=1:points
    for j=1:points
        [SX(i,j),SY(i,j),SZ(i,j)] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u(1),v(j),w(i));
    end
end

surf(SX,SY,SZ,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong');

clear SX SY SZ

 %%%%%%%%%%%%%%%%%%% face 3
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v        u=end  |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->w----
 
SX = zeros(points);
SY = zeros(points);
SZ = zeros(points);

for i=1:points
    for j=1:points
        [SX(i,j),SY(i,j),SZ(i,j)] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u(end),v(j),w(i));
    end
end

surf(SX,SY,SZ,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong');

clear SX SY SZ


 %%%%%%%%%%%%%%%%%%% face 2
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% w               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 
SX = zeros(points);
SY = zeros(points);
SZ = zeros(points);

for i=1:points
    for j=1:points
        [SX(i,j),SY(i,j),SZ(i,j)] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u(i),v(1),w(j));
    end
end

surf(SX,SY,SZ,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong');

clear SX SY SZ

%%%%%%%%%%%%%%%%%%% face 5
%%%%%%%%%%%%%%%%%%%%  ---------------
%%%%%%%%%%%%%%%%%%%% w         v=end |
%%%%%%%%%%%%%%%%%%%% ^               |
%%%%%%%%%%%%%%%%%%%% |               |
%%%%%%%%%%%%%%%%%%%% |               |
%%%%%%%%%%%%%%%%%%%% |               |
%%%%%%%%%%%%%%%%%%%% |               |
%%%%%%%%%%%%%%%%%%%%  --------->u----

SX = zeros(points);
SY = zeros(points);
SZ = zeros(points);

for i=1:points
    for j=1:points
        [SX(i,j),SY(i,j),SZ(i,j)] = VolumePointSpecial(Model,Pwx,Pwy,Pwz,u(i),v(end),w(j));
    end
end

surf(SX,SY,SZ,'edgecolor','none','facecolor',Model.PlotProp.RGB,'FaceLighting','phong');

clear SX SY SZ

if cpoints
    for i = 1 : size(PX,3)
        plot3(PX(:,:,i),PY(:,:,i),PZ(:,:,i),'-','color',Model.PlotProp.ControlRGB,'LineWidth',Model.PlotProp.LineSize); hold on;
        plot3((PX(:,:,i))',(PY(:,:,i))',(PZ(:,:,i))','-','color',Model.PlotProp.ControlRGB,'LineWidth',Model.PlotProp.LineSize);
    end
    
    for i = 1 : size(PX,2)
        PX_(:,:)=PX(:,i,:);
        PY_(:,:)=PY(:,i,:);
        PZ_(:,:)=PZ(:,i,:);
        plot3(PX_,PY_,PZ_,'color',Model.PlotProp.ControlRGB,'LineWidth',Model.PlotProp.LineSize);
        plot3(PX_',PY_',PZ_','color',Model.PlotProp.ControlRGB,'LineWidth',Model.PlotProp.LineSize);
    end
    clear PX_ PY_ PZ_
    for i = 1 : size(PX,1)
        PX_(:,:)=PX(i,:,:);
        PY_(:,:)=PY(i,:,:);
        PZ_(:,:)=PZ(i,:,:);
        plot3(PX_,PY_,PZ_,'color',Model.PlotProp.ControlRGB,'LineWidth',Model.PlotProp.LineSize);
        plot3(PX_',PY_',PZ_','color',Model.PlotProp.ControlRGB,'LineWidth',Model.PlotProp.LineSize);
    end
    
end
 
end
