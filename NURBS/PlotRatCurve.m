function  [ ] = PlotRatCurve(n,p,U,Pw,CPTS,isolines,PlotStruct)


if isolines
    
    Unique = unique(U);
    spans = numel(Unique)-1;
    uu = zeros(spans,floor(1000/spans));
    for i=1:spans
        uu(i,:) = linspace(Unique(i),Unique(i+1),floor(1000/spans));
    end
    C(1:spans,1:floor(1000/spans)) = CPOINT(0,0,0,0,1);
    
    for i=1:spans
        for j=1:numel(uu(1,:))
            C(i,j) = CCurvePoint2(n,p,U,Pw,uu(i,j));
        end
    end
    hold on;
    for i=1:spans
        plot3([C(i,:).x],[C(i,:).y],[C(i,:).z],'LineWidth',PlotStruct.LineSize);
    end
    
    el = unique(U);
    EL(1:numel(el)) = CPOINT(0,0,0,0,1);

    for i=1:numel(el)
        EL(i) = CCurvePoint2(n,p,U,Pw,el(i));
    end
    
    plot3([EL.x],[EL.y],[EL.z],'s','MarkerEdgeColor',PlotStruct.MarkerCPTRGB, 'MarkerFaceColor',PlotStruct.MarkerCPTRGB);
    
else
    
    
    uu = linspace(U(1),U(end),1000);
    C(1:numel(uu)) = CPOINT(0,0,0,0,1);
    
    for i=1:numel(uu)
    	C(i) = CCurvePoint2(n,p,U,Pw,uu(i));
    end
    
    hold on;
    plot3([C.x],[C.y],[C.z],'color',PlotStruct.RGB,'LineWidth',PlotStruct.LineSize);

end

if CPTS == 1 %%wants to plot CPTS
    
    for i=1:numel(Pw)
        wx(i) = Pw(i).x/Pw(i).w;
        wy(i) = Pw(i).y/Pw(i).w;
        wz(i) = Pw(i).z/Pw(i).w;
    end
    
    plot3(wx,wy,wz,'--','color',PlotStruct.ControlRGB);
    plot3(wx,wy,wz,PlotStruct.MarkerCPT,'MarkerEdgeColor',PlotStruct.MarkerCPTRGB, 'MarkerFaceColor',PlotStruct.MarkerCPTRGB);
end


end