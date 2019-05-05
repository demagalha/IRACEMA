function [nq,UQ,mq,VQ,Qw] = SurfaceKnotIns(np,p,UP,mp,q,VP,Pw,dir,uv,k,s,r)

if (dir == 1)
    ump = np+p+1;
    nq = np+r;
    VQ = VP;
    mq = mp;
    Qw(1:np+2,1:mp+1) = CPOINT(0,0,0,0,1);
    
    for i=0:k
        UQ(i+1) = UP(i+1);
    end
    for i=1:r
        UQ(k+i+1) = uv;
    end
    for i=k+1:ump
        UQ(i+r+1) = UP(i+1);
    end
    
    for j=1:r
        L = k-p+j;
        for i=0:p-j+s
            alpha(i+1,j+1) = (uv-UP(L+i+1))/(UP(i+k+2)-UP(L+i+1));
        end
    end
    
    for row=0:mp
        for i=0:k-p
            Qw(i+1,row+1) = Pw(i+1,row+1);
        end
        for i=k-s:np
            Qw(i+r+1,row+1) = Pw(i+1,row+1);
        end
        
        for i=0:p-s
            Rw(i+1) = Pw(k-p+i+1,row+1);
        end
        
        for j=1:r
            L = k-p+j;
            for i=0:p-j-s
                Rw(i+1) = alpha(i+1,j+1)*Rw(i+2) + (1-alpha(i+1,j+1))*Rw(i+1);
            end
            Qw(L+1,row+1) = Rw(1);
            Qw(k+r-j-s+1,row+1) = Rw(p-j-s+1);
        end
        
        for i=L+1:k-s
            Qw(i+1,row+1) = Rw(i-L+1);
        end
    end
end

if(dir == 2)
    vmp = mp+q+1;
    mq = mp+r;
    UQ = UP;
     nq = np;
     Qw(1:np+1,1:mp+2) = CPOINT(0,0,0,0,1);
    
    for i=0:k
        VQ(i+1) = VP(i+1);
    end
    for i=1:r
        VQ(k+i+1) = uv;
    end
    for i=k+1:vmp
        VQ(i+r+1) = VP(i+1);
    end
    
    for j=1:r
        L = k-q+j;
        for i=0:q-j-s
            alpha(i+1,j+1) = (uv-VP(L+i+1))/(VP(i+k+2)-VP(L+i+1));
        end
    end
    
    for row=0:np
        for i=0:k-q
            Qw(row+1,i+1) = Pw(row+1,i+1);
        end
        for i=k-s:mp
            Qw(row+1,i+r+1) = Pw(row+1,i+1);
        end
        
        for i=0:q-s
            Rw(i+1) = Pw(row+1,k-q+i+1);
        end
        
        for j=1:r
            L = k-q+j;
            for i=0:q-j-s
                Rw(i+1) = alpha(i+1,j+1)*Rw(i+2) + (1-alpha(i+1,j+1))*Rw(i+1);
            end
            Qw(row+1,L+1) = Rw(1);
            Qw(row+1,k+r-j-s+1) = Rw(q-j-s+1);
        end
        
        for i=L+1:k-s
            Qw(row+1,i+1) = Rw(i-L+1);
        end
    end
end
end