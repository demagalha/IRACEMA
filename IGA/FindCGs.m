% function [tau, tav, taw] = FindCGs(Model)

U = Model.U;
pu = Model.pu;
V = Model.V;
pv = Model.pv;
W = Model.W;
pw = Model.pw;
p = [pu, pv, pw];
Knots = {U, V, W};
dN = {};
initial ={};
CG = {};
% Initial Guess for the Collocation Points
greville = cell(3,1);
tmp = 0;
for i=1:3
        for j=1:p(i)+2
            for k=j+1:j+p(i)
%             greville{i}(j) = greville{i}(j) + Knots{i}(j)/p(i);
               tmp = tmp +Knots{i}(k)/p(i);
            end
            greville{i}(j) = tmp;
            tmp = 0;
        end
        greville{i} = unique(greville{i});
end
% for i=1:3
%     initial = greville{i};
%     for j=1:numel(initial)
%         s = FindSpanLinear(length(Knots{i})-p(i)-1,p(i),initial(j),Knots{i});
%         dN = DersBasisFun(s,initial(j),p(i),3,U);
%     tol = 1;
%     ctr = 0;
%         while (tol>1e-3 && ctr < 5)
%            tol = initial(j);
%            initial(j) = initial(j) - sum(dN(1,:).*dN(3,:))/sum(dN(1,:).*dN(4,:)); 
%            tol = abs(tol - initial(j));
%            ctr = ctr+1
%         end
%     end
%     CG{i} = initial;
% end


