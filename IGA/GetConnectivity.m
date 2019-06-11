function [INN, IEN, nel, nen] = GetConnectivity(varargin)
    switch nargin
        case 6
            nu = varargin{1};
            pu = varargin{2};
            nv = varargin{3};
            pv = varargin{4};
            nw = varargin{5};
            pw = varargin{6};
            
            nu = nu+1;
            nv = nv+1;
            nw = nw+1;

            nel = (nu-pu)*(nv-pv)*(nw-pw); %number of elements
            nnp = nu*nv*nw;   %%number of global basis funs
            nen = (pu+1)*(pv+1)*(pw+1); %%number of local basis funs
            INN = zeros(nnp,3);
            IEN = zeros(nen,nel);

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
        case 4
            nu = varargin{1};
            pu = varargin{2};
            nv = varargin{3};
            pv = varargin{4};
            nu = nu+1;
            nv = nv+1;
            nel = (nu-pu)*(nv-pv); %number of elements
            nnp = nu*nv;   %%number of global basis funs
            nen = (pu+1)*(pv+1); %%number of local basis funs
            INN=zeros(nnp,2);
            IEN=zeros(nen,nel);

            A = 0;
            e = 0;

               for j=1:nv
                    for i=1:nu
                        A = A +1;

                        INN(A,1) = i;
                        INN(A,2) = j;
                        if i >= pu+1 && j >= pv+1
                            e = e+1;

                                 for jloc=0:pv
                                    for iloc=0:pu
                                        B = A - jloc*nu -iloc;

                                        b = jloc*(pu+1) + iloc +1;

                                        IEN(b,e) = B;
                                    end
                                end
                            end
                        end
               end
        case 2
            nu = varargin{1}+1;
            pu = varargin{2};
            nel = (nu-pu); %number of elements
            nnp = nu;   %%number of global basis funs
            nen = (pu+1); %%number of local basis funs
            INN =zeros(nnp,1);
            IEN=zeros(nen,nel);
            A = 0;
            e = 0;

       for i=1:nu
           A = A +1;
           INN(A,1) = i;
           if i >= pu+1
              e = e+1;
              for iloc=0:pu
                  B = A -iloc;
                  b = iloc +1;
                  IEN(b,e) = B;
              end
          end
       end
    end
end