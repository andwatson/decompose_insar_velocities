function gnss_grd2mat(in_north,in_east,in_up,outfile,error_north,error_east,error_up)
%% gnss_grd2mat.m
% Loads three .grd file contains North, East and vertical velocities (made with gmt)
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
[~,~,up] = grdread2(in_up);

%% format

gnss_field.x = x;
gnss_field.y = y;
gnss_field.N = north;
gnss_field.E = east;
gnss_field.U = up;

%% optionally load uncertainties

if nargin == 5
    
    [~,~,sN] = grdread2(error_north);
    [~,~,sE] = grdread2(error_east);
    [~,~,sU] = grdread2(error_up);
    
    gnss_field.sN = sN;
    gnss_field.sE = sE;
    gnss_field.sU = sU;
    
end

%% output

save(outfile,'gnss_field')

end
