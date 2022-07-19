%% compare_to_gnss.m
%
% Compare a veloicity field (stored as a tif) to GNSS station velocities.
%
% Andrew Watson     06-06-2022

%% setup

% input files
vel_dir = '/scratch/eearw/decomp_frame_vels/out/ref_tests/';
vel_file = 'filt_81km_vE.geo.tif';
gnss_file = '/nfs/a285/homes/eearw/gmt/khorrami_gnss_trimmed.csv';

% component toggle - options are East, North, LOS
comp = 'East';

%% load inputs

% load geotifs
[lon,lat,vel,~,~] = read_geotiff([vel_dir vel_file]);

% load gnss vels
gnss = readmatrix(gnss_file);

%% calculate residuals

% isolate gnss component for comparison
switch comp
    case 'East'
        gnss_comp = gnss(:,[1 2 3]);
        
    case 'North'
        gnss_comp = gnss(:,[1 2 4]);
        
end

% calculate residual
gnss_resid = zeros(size(gnss_comp,1),1);

for ii = 1:size(gnss_comp,1)
    
    % index of closest vel to gnss stations
    [~,ind_x] = min(abs(lon-gnss_comp(ii,1)));
    [~,ind_y] = min(abs(lat-gnss_comp(ii,2)));
    
    gnss_resid(ii) = gnss_comp(ii,3) - vel(ind_y,ind_x);
end

%% plot results

% set plotting parameters
lonlim = [min(lon) max(lon)];
latlim = [min(lat) max(lat)];
clim = [-10 10];
load('/nfs/a285/homes/eearw/gmt/colourmaps/vik/vik.mat')

figure(); hold on
plt_data(lon,lat,vel,lonlim,latlim,clim,'East and GNSS residuals (mm/yr)',[],[])
scatter(gnss_comp(:,1),gnss_comp(:,2),70,gnss_resid,'Filled','MarkerEdgeColor','k')
colormap(vik)

figure()
histogram(gnss_resid,20);
title('GNSS residual (mm/yr)')










% [~,~,compE,~,~] = read_geotiff([vel_dir Efile]);
% [~,~,compN,~,~] = read_geotiff([vel_dir Nfile]);
% [~,~,compU,~,~] = read_geotiff([vel_dir Ufile]);

% case 'LOS'
%     gnss_comp = gnss(:,1:3);
% 
%     for ii = 1:size(gnss_comp,1)
% 
%         % get index
%         [~,ind_x] = min(abs(lon-gnss_los(ii,1)));
%         [~,ind_y] = min(abs(lat-gnss_los(ii,2)));
% 
%         % project
%         gnss_comp(ii,3) = gnss(ii,3).*compE(ind_y,ind_x) ...
%             + gnss(ii,4).*compN(ind_y,ind_x);
% 
%     end