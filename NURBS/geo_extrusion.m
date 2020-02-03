function extrusion = geo_extrusion(Model,vector)

%Model to be extruded
%vector with the direction to be extruded. e.g. [0,0,20], which would
%extrude Model by 20 units in the z direction

S = Model.get_point_cell;
U = Model.U;
V = Model.V;
nu = Model.nu;
pu = Model.pu;
nv = Model.nv;
pv = Model.pv;

vec = [vector 0];
B = S;
for i = 1:nu+1
    for j = 1:nv+1
        B{i,j,2} = B{i,j,1} + vec;
    end
end

extrusion = Geometry('volume',pu,U,pv,V,1,[0,0,1,1],B);

end