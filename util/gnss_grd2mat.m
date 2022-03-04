function gnss_grd2mat(in_north,in_east,outfile)
%% gnss_grd2mat.m
% Loads two .grd file contains North and East velocities (made with gmt)
% and outputs a mat file that can be used in decompose_insar_velocities.
% Inputs must share the same grid.
%
% Andrew Watson     2022-02-16

%% load

[x,y,north] = grdread2(in_north);
[~,~,east] = grdread2(in_east);

%% format and output

gnss_field.x = x;
gnss_field.y = y;
gnss_field.N = north;
gnss_field.E = east;

save(outfile,'gnss_field')

end