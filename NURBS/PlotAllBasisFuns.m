function [ ] = PlotAllBasisFuns(n,p,U)

%%%%plot all basis functions in a single direction
%%input n (n+1 ctrl points), p (degree),U (knot vector)

uu = linspace(0,U(end),1000);
m = n+p+1;
Nip = zeros(n+1,numel(uu));
for i=1:numel(uu)
    for j=1:n+1
        Nip(j,i) = OneBasisFun(p,m,U,j-1,uu(i));
    end
end

plot(uu,Nip(1,:));
hold on;
for i=2:n+1
    plot(uu,Nip(i,:));
end

end