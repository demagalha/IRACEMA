function [ ] = PlotOneBasisFuns(n,p,U,j)

%%%%plot a single basis functions in a single direction
%%input n (n+1 ctrl points), p (degree),U (knot vector)

uu = linspace(0,U(end),1000);
m = n+p+1;
Nip = zeros(1,numel(uu));
for i=1:numel(uu)
        Nip(1,i) = OneBasisFun(p,m,U,j-1,uu(i));
end

Nip(Nip == 0) = nan;
plot(uu,Nip(1,:));
hold on;

end