function [] = PlotRatVolElement(Model,span1,span2,span3)

Unique = unique(Model.U);
Vnique = unique(Model.V);
Wnique = unique(Model.W);

u = linspace(Unique(span1),Unique(span1+1),10);
v = linspace(Vnique(span2),Vnique(span2+1),10);
w = linspace(Wnique(span3),Wnique(span3+1),10);

 %%%%%%%%%%%%%%%%%%% face 1
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 
 SX = zeros(numel(u),numel(v));
 SY = zeros(numel(u),numel(v));
 SZ = zeros(numel(u),numel(v));
 
 for i=1:numel(u)
    for j=1:numel(v)
        point = Model.eval_point(u(i),v(j),w(1));
        SX(i,j) = point.x;
        SY(i,j) = point.y;
        SZ(i,j) = point.z;
        %point.z
    end
 end
 
 h = surf(SX,SY,SZ);
 set(h,'edgecolor','none','facecolor',[1 0 0])
 
 hold on;
 clear SX SY SZ point
 
 %%%%%%%%%%%%%%%%%%% face 6
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v         w=end |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 

 SX = zeros(numel(u),numel(v));
 SY = zeros(numel(u),numel(v));
 SZ = zeros(numel(u),numel(v));
 
 for i=1:numel(u)
    for j=1:numel(v)
        point = Model.eval_point(u(i),v(j),w(end));
        SX(i,j) = point.x;
        SY(i,j) = point.y;
        SZ(i,j) = point.z;
        %point.z
    end
 end
 
 h = surf(SX,SY,SZ);
 set(h,'edgecolor','none','facecolor',[0 1 0])
 clear SX SY SZ point
 
 %%%%%%%%%%%%%%%%%%% face 4
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->w----
 
 SX = zeros(numel(w),numel(v));
 SY = zeros(numel(w),numel(v));
 SZ = zeros(numel(w),numel(v));
 
 for i=1:numel(w)
    for j=1:numel(v)
        point = Model.eval_point(u(1),v(j),w(i));
        SX(i,j) = point.x;
        SY(i,j) = point.y;
        SZ(i,j) = point.z;
    end
 end
 
 h = surf(SX,SY,SZ);
 set(h,'edgecolor','none','facecolor',[0 0 1])
 clear SX SY SZ point
 
 
 %%%%%%%%%%%%%%%%%%% face 3
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% v        u=end  |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->w----
 
 SX = zeros(numel(w),numel(v));
 SY = zeros(numel(w),numel(v));
 SZ = zeros(numel(w),numel(v));
 
 for i=1:numel(w)
    for j=1:numel(v)
        point = Model.eval_point(u(end),v(j),w(i));
        SX(i,j) = point.x;
        SY(i,j) = point.y;
        SZ(i,j) = point.z;
    end
 end
 
 h = surf(SX,SY,SZ);
 set(h,'edgecolor','none','facecolor',[1 1 0])
 clear SX SY SZ point

 %%%%%%%%%%%%%%%%%%% face 2
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% w               |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 
 SX = zeros(numel(u),numel(w));
 SY = zeros(numel(u),numel(w));
 SZ = zeros(numel(u),numel(w));
 
 for i=1:numel(u)
    for j=1:numel(w)
        point = Model.eval_point(u(i),v(1),w(j));
        SX(i,j) = point.x;
        SY(i,j) = point.y;
        SZ(i,j) = point.z;
    end
 end
 
 h = surf(SX,SY,SZ);
 set(h,'edgecolor','none','facecolor',[0 1 1])
 clear SX SY SZ point
 
 %%%%%%%%%%%%%%%%%%% face 5
 %%%%%%%%%%%%%%%%%%%%  ---------------
 %%%%%%%%%%%%%%%%%%%% w         v=end |
 %%%%%%%%%%%%%%%%%%%% ^               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%% |               |
 %%%%%%%%%%%%%%%%%%%%  --------->u----
 
 SX = zeros(numel(u),numel(w));
 SY = zeros(numel(u),numel(w));
 SZ = zeros(numel(u),numel(w));
 
 for i=1:numel(u)
    for j=1:numel(w)
        point = Model.eval_point(u(i),v(end),w(j));
        SX(i,j) = point.x;
        SY(i,j) = point.y;
        SZ(i,j) = point.z;
    end
 end
 
 h = surf(SX,SY,SZ);
 set(h,'edgecolor','none','facecolor',[0.2 0.2 0.2])
 clear SX SY SZ point
 
 
end
 

