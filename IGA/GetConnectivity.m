function [INN, IEN, nel, nen] = GetConnectivity(nu,pu,nv,pv,nw,pw)

nu = nu+1;
nv = nv+1;
nw = nw+1;

nel = (nu-pu)*(nv-pv)*(nw-pw); %number of elements
nnp = nu*nv*nw;   %%number of global basis funs
nen = (pu+1)*(pv+1)*(pw+1); %%number of local basis funs
INN(1:nnp,1:3) = 0;
IEN(1:nen,1:nel) = 0;

A = 0;
e = 0;

for k=1:nw
    for j=1:nv
        for i=1:nu
            A = A +1;
            
            INN(A,1) = i;
            INN(A,2) = j;
            INN(A,3) = k;
            
            if i >= pu+1 && j >= pv+1 && k >= pw+1
                e = e+1;
                
                for kloc=0:pw
                    for jloc=0:pv
                        for iloc=0:pu
                            B = A - kloc*nu*nv - jloc*nu -iloc;
                            
                            b = kloc*(pu+1)*(pv+1) + jloc*(pu+1) + iloc +1;
                            
                            IEN(b,e) = B;
                        end
                    end
                end
            end
        end
    end
end
end