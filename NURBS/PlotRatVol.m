function [] = PlotRatVol(PX,PY,PZ,w,nu,pu,U,nv,pv,V,nw,pw,W,render)

Unique = unique(U);
Vnique = unique(V);
Wnique = unique(W);

 %%%%%%%%%%%%%%%%%%%% plotar as faces
 %%%%%%%%%%%%%%%%%%% face 1
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 

FX(:,:) = PX(:,:,1);
FY(:,:) = PY(:,:,1);
FZ(:,:) = PZ(:,:,1);
FW(:,:) = w(:,:,1);

PlotSurf(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render);
    hold on;
for i=1:numel(Unique)
    PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Unique(i));
end

for j=1:numel(Vnique)
    PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Vnique(j));
end
        
 %%%%%%%%%%%%%%%%%%% face 6
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v         w=end |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----

FX(:,:) = PX(:,:,end);
FY(:,:) = PY(:,:,end);
FZ(:,:) = PZ(:,:,end);
FW(:,:) = w(:,:,end);

PlotSurf(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render);
for i=1:numel(Unique)
    PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Unique(i));
end

for j=1:numel(Vnique)
   PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Vnique(j));
end

 %%%%%%%%%%%%%%%%%%% face 4
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->w----
 
 clear FX FY FZ FW
 FX(:,:) = PX(1,:,:);
 FX = FX';
 FY(:,:) = PY(1,:,:);
 FY = FY';
 FZ(:,:) = PZ(1,:,:);
 FZ = FZ';
 FW(:,:) = w(1,:,:);
 FW = FW';
 
 
 PlotSurf(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render);
 for i=1:numel(Wnique)
    PlotSurfPatches_1(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Wnique(i));
end

for j=1:numel(Vnique)
    PlotSurfPatches_2(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Vnique(j));
end
     
 %%%%%%%%%%%%%%%%%%% face 3
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v        u=end  |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->w----
 
 clear FX FY FZ FW
 FX(:,:) = PX(end,:,:);
 FX = FX';
 FY(:,:) = PY(end,:,:);
 FY = FY';
 FZ(:,:) = PZ(end,:,:);
 FZ = FZ';
 FW(:,:) = w(end,:,:);
 FW = FW';
 
PlotSurf(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render);
for i=1:numel(Wnique)
    PlotSurfPatches_1(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Wnique(i));
end

for j=1:numel(Vnique)
    PlotSurfPatches_2(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Vnique(j));
end
     
 %%%%%%%%%%%%%%%%%%% face 2
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% w               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----

 clear FX FY FZ FW
 FX(:,:) = PX(:,1,:);
 FY(:,:) = PY(:,1,:);
 FZ(:,:) = PZ(:,1,:);
 FW(:,:) = w(:,1,:);
 
 PlotSurf(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render);
 for i=1:numel(Unique)
    PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Unique(i));
end

for j=1:numel(Wnique)
    PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Wnique(j));
end
     
 %%%%%%%%%%%%%%%%%%% face 5
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% w         v=end |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
    
 clear FX FY FZ FW
 FX(:,:) = PX(:,end,:);
 FY(:,:) = PY(:,end,:);
 FZ(:,:) = PZ(:,end,:);
 FW(:,:) = w(:,end,:);
 
 PlotSurf(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render);
for i=1:numel(Unique)
    PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Unique(i));
end

for j=1:numel(Wnique)
    PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Wnique(j));
end
 light;

 
for i = 1 : size(PX,3)
plot3(PX(:,:,i),PY(:,:,i),PZ(:,:,i),'-','color','black','LineWidth',2); hold on;
plot3((PX(:,:,i))',(PY(:,:,i))',(PZ(:,:,i))','-','color','black','LineWidth',2);
end
for i = 1 : size(PX,2)
PX_(:,:)=PX(:,i,:);
PY_(:,:)=PY(:,i,:);
PZ_(:,:)=PZ(:,i,:);
plot3(PX_,PY_,PZ_,'color','black','LineWidth',2);
plot3(PX_',PY_',PZ_','color','black','LineWidth',2);
end
clear PX_ PY_ PZ_
for i = 1 : size(PX,1)
PX_(:,:)=PX(i,:,:);
PY_(:,:)=PY(i,:,:);
PZ_(:,:)=PZ(i,:,:);
plot3(PX_,PY_,PZ_,'color','black','LineWidth',2);
plot3(PX_',PY_',PZ_','color','black','LineWidth',2);
end
 
end 
