function [nnew, Unew, Qw] = DegreeElevateSolid(nu,pu,U,nv,pv,V,nw,pw,W,Pw,dir,t)

%[nh,Uh,Qw] = DegreeElevateCurve(n,p,U,Pw,t)
%[nh,Uh,mh,Vh,Qw] = DegreeElevateSurface(n,p,U,m,q,V,Pw,dir,t)

if dir == 1
    
s = FindS(U,pu);
Qw(1:size(Pw,1)+t*(1+s),1:size(Pw,2),1:size(Pw,3)) = CPOINT(0,0,0,0,0);

P(1:size(Pw,1)) = CPOINT(0,0,0,0,1);

    
    for j=1:size(Pw,2)
        for k=1:size(Pw,3)
            for i=1:size(Pw,1)
                P(i) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            
            [nh,Ubar,Q] = DegreeElevateCurve(nu,pu,U,P,t);
         
            for i=1:size(Pw,1)+t*(1+s)
             Qw(i,j,k) = Q(i);
            end
        end
    end

Unew = Ubar;
nnew = nh;
end

if dir == 2
    
    s = FindS(V,pv);
    Qw(1:size(Pw,1),1:size(Pw,2)+t*(1+s),1:size(Pw,3)) = CPOINT(0,0,0,0,0);
    
    P(1:size(Pw,2)) = CPOINT(0,0,0,0,1);
    
    for i=1:size(Pw,1)
        for k=1:size(Pw,3)
            for j=1:size(Pw,2)
                P(j) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            
            [nh,Vbar,Q] = DegreeElevateCurve(nv,pv,V,P,t);
         
            for j=1:size(Pw,2)+t*(1+s)
             Qw(i,j,k) = Q(j);
            end
        end
    end
    
Unew = Vbar;
nnew = nh;
end

if dir == 3
    
    s = FindS(W,pw);
    Qw(1:size(Pw,1),1:size(Pw,2),1:size(Pw,3)+t*(1+s)) = CPOINT(0,0,0,0,0);
    
    P(1:size(Pw,3)) = CPOINT(0,0,0,0,1);
    
    for i=1:size(Pw,1)
        for j=1:size(Pw,2)
            for k=1:size(Pw,3)
                P(k) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            
            [nh,Wbar,Q] = DegreeElevateCurve(nw,pw,W,P,t);
            %nh
            %size(Pw,3)+t*(1+s)
            %nh+t*(1+s)
         
            for k=1:size(Pw,3)+t*(1+s)
             Qw(i,j,k) = Q(k);
            end
        end
    end
    
Unew = Wbar;
nnew = nh;
end


end
