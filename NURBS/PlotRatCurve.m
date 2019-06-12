function  [ ] = PlotRatCurve(n,p,U,Pw,CPTS,NOZ)


uu = linspace(U(1),U(end),1000);
C(1:numel(uu)) = CPOINT(0,0,0,0,1);

for i=1:numel(uu)
    C(i) = CCurvePoint2(n,p,U,Pw,uu(i));
end

if CPTS == 1 %%wants to plot CPTS
    
    for i=1:numel(Pw)
        wx(i) = Pw(i).x/Pw(i).w;
        wy(i) = Pw(i).y/Pw(i).w;
        wz(i) = Pw(i).z/Pw(i).w;
    end

    if NOZ == 1
        plot([C.x],[C.y],'color','blue');
    else
        plot3([C.x],[C.y],[C.z],'color','blue');
    end
    
    hold on;
    plot3(wx,wy,wz,'--','color','black');
    plot3(wx,wy,wz,'o','MarkerEdgeColor',[1 0 0], 'MarkerFaceColor',[1 0 0]);

else
    if NOZ == 1
        plot([C.x],[C.y],'color','blue');
    else
        plot3([C.x],[C.y],[C.z],'color','blue');
    end
end

end