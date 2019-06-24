clear all
close all
clc
load('plate_modal.mat')
 B = Model.get_point_cell;
 u = cell(size(B));
 comb = u;
 [sz1 sz2] = size(autovector);
Modos = cell(sz2,1);
 for n_modo = 1:sz2
    for i =1:size(ID,2)
        u{i} = [autovector(3*(i-1) +1,n_modo),autovector(3*(i-1) +2,n_modo),autovector(3*(i-1) +3,n_modo), 0];
        comb{i} = B{i} + u{i};
    end
    Modos{n_modo} = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
 end
 
 Modos{1}.plot_geo
