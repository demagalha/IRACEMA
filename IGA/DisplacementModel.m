function D = DisplacementModel(autovector,ID,Model)

P = Model.get_point_cell;
[size_u, size_v, size_w] = size(P);
B = cell(size_u,size_v,size_w);


size_ID = size(ID,2);
for i=1:size_ID
    B{i} = [autovector(ID(1,i)), autovector(ID(2,i)), autovector(ID(3,i)), P{i}(4)];
end

D = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,B);
    


end