function gnss_grd2mat(in_north,in_east,outfile,error_north,error_east)
%% gnss_grd2mat.m
% Loads two .grd file contains North and East velocities (made with gmt)
% and outputs a mat file that can be used in decompose_insar_velocities.
% Inputs must share the same grid.
%
% Optional - loads and saves uncertainties associated with velocities.
%
% Andrew Watson     2022-02-16

%% load

% load velocities
[x,y,north] = grdread2(in_north);
[~,~,east] = grdread2(in_east);

%% format

gnss_field.x = x;
gnss_field.y = y;
gnss_field.N = north;
gnss_field.E = east;

%% optionally load uncertainties

if nargin == 5
    
    [~,~,sN] = grdread2(error_north);
    [~,~,sE] = grdread2(error_east);
    
    gnss_field.sN = sN;
    gnss_field.sE = sE;
    
end

%% output

save(outfile,'gnss_field')

end