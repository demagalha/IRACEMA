% function Rules = getWQ(pu,U)
%% Generate the Quadrature Point Grid:
    Omega = unique(U);
    midpoints = 0.5*(Omega(1:end-1) +Omega(2:end))';
    first = linspace(Omega(1),Omega(2),pu+1);
    last = linspace(Omega(end-1),Omega(end),pu+1);
    qpoints = union(first,last);    
    qpoints = union(qpoints,midpoints);
    qpoints = union(qpoints,Omega);

%% Generate K, I and Q sets
X = unique(U); % Distinct Knots of U.
ndof = length(U)-pu-1;
nel = length(X)-1;
nq = length(qpoints);
K = cell(ndof,1);
I = cell(ndof,1);
Q = cell(ndof,1);
for i=1:ndof
    support = i-1:i+pu-1;
    si = FindSpanLinear(ndof-1,pu,U(i),U);
    for k=1:nel
        span = FindSpanLinear(ndof-1,pu,X(k),U);
        [log, ~] = ismember(support,span);
        if any(log(1:end-1))
            if k==1
                K{i}(1) = k;
            else
                K{i} = [K{i} k];
            end
        end
    end
    for j=1:ndof
        troppus = j-1:j+pu-1;
        overlap = intersect(support,troppus);
        if ~isempty(overlap)
            if j == 1
                I{i}(1) = j;
            else
                I{i} = [I{i} j];
            end
        end
    end
    for q=1:nq
        span = FindSpanLinear(ndof-1,pu,qpoints(q),U);
        [log, ~] = ismember(support,span);
        if any(log)
            if q==1
                Q{i}(1) = q;
            else
                Q{i} = [Q{i} q];
            end
        end
    end
end
        
% end