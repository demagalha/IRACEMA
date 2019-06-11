a = mfilename('fullpath');
IGA = strrep(a,'startup','IGA');
Examples = strrep(a,'startup','Examples');
NURBS = strrep(a,'startup','NURBS');

addpath (IGA)
addpath (Examples)
addpath (NURBS)

clear a IGA Examples NURBS

%ideally we would want to add to the matlab search path \MATLAB\...\bin a
%simple startup file or add the whole package inside the searchpath, this
%file is just a work around for it, since IRACEMA is far from being finished