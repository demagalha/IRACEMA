function Rules = getWQ(p, knots)
    sz = length(p);
    Rules = cell(sz,1);
    for d=1:sz
    ndof = length(knots{d})-p(d)-1;
    distinct_knots = unique(knots{d});
    nel = numel(distinct_knots) -1; % Page 7, Calabrò
    q_first_el = linspace(distinct_knots(1), distinct_knots(2), p(d)+2);
    q_last_el = linspace(distinct_knots(end-1), distinct_knots(end), p(d)+2);
    q_int_el = sort([distinct_knots(3:end-2) 0.5*(distinct_knots(2:end-2)+distinct_knots(3:end-1))]);
    all_points = [q_first_el q_int_el q_last_el];
    all_points = unique(all_points);
    nq = length(all_points);
    % Construction of sets K, I and Q
    % K is composed by integers k such that the span ]X_k-1, X_k[ is inside
    % the support of B_il.
    % I is composed of the integers j such that B_i (dot) B_j ~= 0
    % Q is composed of integers q such that x_q lies in the supp of B_i
    K = cell(1,ndof);
    I = cell(1,ndof);
    Q = cell(1,ndof);
    for i=1:ndof
      for k=1:nel
            s = FindSpanLinear(ndof-1,p(d),distinct_knots(k),knots{d});
            sup = i-1:i+p(d)-1;
            [log, ~] = ismember(sup,s);
            if any(log)
                K{i} = [K{i} k];
            end
      end
      K{i} = unique(K{i});
      for j=1:ndof
          sup = i:i+p(d);
          sup2 = j:j+p(d);
          dummy = intersect(sup,sup2);
          if isempty(dummy)
              continue
          else
              I{i} = [I{i} j];
          end
      end
      I{i} = unique(I{i});
      for q=1:nq
          sup = i-1:i+p(d)-1;
          s = FindSpanLinear(ndof-1,p(d),all_points(q),knots{d});
          [log, ~] = ismember(sup,s);
          if any(log)
              Q{i} = [Q{i} q];
          end
      end
      Q{i} = unique(Q{i});
    end
 %% Algorithm 1 of the paper, Initializations
    B0 = zeros(ndof,nq);
    B1 = B0;
    for i=1:ndof
        for q=1:nq
            s = FindSpanLinear(ndof-1,p(d),all_points(q),knots{d});
            supp = i-1:i+p(d)-1;
            [log, ~] = ismember(supp,s);
            if any(log)
                N = DersBasisFun(s,all_points(q),p(d),1,knots{d});
                B0(i,q) = N(1,end-(s+1)+i);
                B1(i,q) = N(2,end-(s+1)+i);
            end
        end
    end
    Rules{d}.Basis = B0;
    Rules{d}.DerBasis = B1;
    % Integrals
    [x_gauss, w_gauss] = getGP(p(d)+1);
    % from [-1,1] to [0,1];
    x_gauss = x_gauss/2 +1/2;
    w_gauss = w_gauss/2;
    % Normalizing the knot vector
    nknots = knots{d}/knots{d}(end);
    I00 = zeros(ndof,ndof);
    I01 = zeros(ndof,ndof);
    I10 = zeros(ndof,ndof);
    I11 = zeros(ndof,ndof);
    for i=1:ndof
        Interval = I{i};
        for j=Interval(1:end)
            for q=1:length(x_gauss)
                supp1 = i-1:i+p(d)-1;
                supp2 = j-1:j+p(d)-1;
                s = FindSpanLinear(ndof-1,p(d),x_gauss(q),knots{d});
                [log1, ~] = ismember(supp1,s);
                [log2, ~] = ismember(supp2,s);
                if any(log1)
                    N = DersBasisFun(s,x_gauss(q),p(d),1,knots{d});
                    Bi = N(1,end-(s+1)+i);
                    dBi = N(2,end-(s+1)+i);
                else
                    Bi = 0;
                    dBi = 0;
                end
                if any(log2)
                    N = DersBasisFun(s,x_gauss(q),p(d),1,knots{d});
                    Bj = N(1,end-(s+1)+j);
                    dBj = N(2,end-(s+1)+j);
                else
                    Bj = 0;
                    dBj = 0;
                end
                I00(i,j) = I00(i,j) + Bi*Bj*w_gauss(q);
                I01(i,j) = I01(i,j) +Bi*dBj*w_gauss(q);
                I10(i,j) = I10(i,j) +dBi*Bj*w_gauss(q);
                I11(i,j) = I11(i,j) +dBi*dBj*w_gauss(q);
            end
        end
    end

%% Algorithm 2 - Construction of Univariate WQ rules
w00 = cell(ndof,1);
w01 = w00;
w10 = w00;
w11 = w00;
for i=1:ndof
   w00{i} = lsqminnorm(B0(I{i},Q{i}),I00(i,I{i})');
   w10{i} = lsqminnorm(B1(I{i},Q{i}),I10(i,I{i})');
   w01{i} = lsqminnorm(B0(I{i},Q{i}),I01(i,I{i})');
   w11{i} = lsqminnorm(B1(I{i},Q{i}),I11(i,I{i})');
end
Weights = cell(4,1);
Weights{1} = w00;
Weights{2} = w01;
Weights{3} = w10;
Weights{4} = w11;
Rules{d}.Weights = Weights;
Rules{d}.I = I;
Rules{d}.Q = Q;
Rules{d}.K = K;
Rules{d}.Points = all_points;
Rules{d}.ndof = ndof;
end
end