a = mfilename('fullpath');
IGA = strrep(a,'startup','IGA');
Examples = strrep(a,'startup','Examples');
NURBS = strrep(a,'startup','NURBS');
Models = strrep(a,'startup','Models');
UserGuide = strrep(a,'startup','UserGuide');

addpath (IGA)
addpath (Examples)
addpath (NURBS)
addpath (Models)
addpath (UserGuide)
clear a IGA Examples NURBS Models UserGuide

%ideally we would want to add to the matlab search path \MATLAB\...\bin a
%simple startup file or add the whole package inside the searchpath, this
%file is just a work around for it, since IRACEMA is far from being finished