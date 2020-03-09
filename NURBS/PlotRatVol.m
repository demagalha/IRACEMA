function [] = PlotRatVol(Model,render,cpoints,isolines)

Unique = unique(Model.U);
Vnique = unique(Model.V);
Wnique = unique(Model.W);

PX = Model.PX; PY = Model.PY; PZ = Model.PZ; w = Model.weight;
nu = Model.nu; nv = Model.nv; nw = Model.nw;
pu = Model.pu; pv = Model.pv; pw = Model.pw;
U = Model.U; V = Model.V; W = Model.W;
PlotStruct = Model.PlotProp;

P = ExtractFaces(Model);
F1 = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,P{1}.P);
F6 = Geometry('surf',Model.pu,Model.U,Model.pv,Model.V,P{6}.P);

F4 = Geometry('surf',Model.pw,Model.W,Model.pv,Model.V,P{4}.P);
F3 = Geometry('surf',Model.pw,Model.W,Model.pv,Model.V,P{3}.P);

F2 = Geometry('surf',Model.pu,Model.U,Model.pw,Model.W,P{2}.P);
F5 = Geometry('surf',Model.pu,Model.U,Model.pw,Model.W,P{5}.P);

%
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

PlotSurf(F1,render);
    hold on;
    
if isolines
    for i=1:numel(Unique)
            PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Unique(i),PlotStruct);
    end
    
    for j=1:numel(Vnique)
        PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Vnique(j),PlotStruct);
    end
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

PlotSurf(F6,render);

if isolines
    for i=1:numel(Unique)
        PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Unique(i),PlotStruct);
    end

    for j=1:numel(Vnique)
        PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nv,pv,V,render,Vnique(j),PlotStruct);
    end
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
 
 
PlotSurf(F4,render);

if isolines
    
    for i=1:numel(Wnique)
        PlotSurfPatches_1(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Wnique(i),PlotStruct);
    end

    for j=1:numel(Vnique)
        PlotSurfPatches_2(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Vnique(j),PlotStruct);
    end
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
 
PlotSurf(F3,render);

if isolines
    for i=1:numel(Wnique)
        PlotSurfPatches_1(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Wnique(i),PlotStruct);
    end

    for j=1:numel(Vnique)
        PlotSurfPatches_2(FX,FY,FZ,FW,nw,pw,W,nv,pv,V,render,Vnique(j),PlotStruct);
    end
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
 
PlotSurf(F2,render);

if isolines
    for i=1:numel(Unique)
        PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Unique(i),PlotStruct);
    end

    for j=1:numel(Wnique)
        PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Wnique(j),PlotStruct);
    end
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
 
PlotSurf(F5,render);

if isolines    
    for i=1:numel(Unique)
        PlotSurfPatches_1(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Unique(i),PlotStruct);
    end

    for j=1:numel(Wnique)
        PlotSurfPatches_2(FX,FY,FZ,FW,nu,pu,U,nw,pw,W,render,Wnique(j),PlotStruct);
    end
end

light;
%
if cpoints
    for i = 1 : size(PX,3)
        plot3(PX(:,:,i),PY(:,:,i),PZ(:,:,i),'-','color',PlotStruct.ControlRGB,'LineWidth',PlotStruct.LineSize); hold on;
        plot3((PX(:,:,i))',(PY(:,:,i))',(PZ(:,:,i))','-','color',PlotStruct.ControlRGB,'LineWidth',PlotStruct.LineSize);
    end
    
    for i = 1 : size(PX,2)
        PX_(:,:)=PX(:,i,:);
        PY_(:,:)=PY(:,i,:);
        PZ_(:,:)=PZ(:,i,:);
        plot3(PX_,PY_,PZ_,'color',PlotStruct.ControlRGB,'LineWidth',PlotStruct.LineSize);
        plot3(PX_',PY_',PZ_','color',PlotStruct.ControlRGB,'LineWidth',PlotStruct.LineSize);
    end
    clear PX_ PY_ PZ_
    for i = 1 : size(PX,1)
        PX_(:,:)=PX(i,:,:);
        PY_(:,:)=PY(i,:,:);
        PZ_(:,:)=PZ(i,:,:);
        plot3(PX_,PY_,PZ_,'color',PlotStruct.ControlRGB,'LineWidth',PlotStruct.LineSize);
        plot3(PX_',PY_',PZ_','color',PlotStruct.ControlRGB,'LineWidth',PlotStruct.LineSize);
    end
    
end
 
end 
