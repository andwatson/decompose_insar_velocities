function [vel] = plate_motion_bias(par,x,y,vel,compE,compN)
%=================================================================
% function plate_motion_bias()
%-----------------------------------------------------------------
% Mitigate the "reference frame effect" caused by rigid plate motions by
% subtracting the No-Net-Rotation plate motion in ITRF from the LOS
% velocities.
%
% The input plate motion vector file should contain (at least) the
% following columns in positions 1-4 : lon, lat, East vel, North vel.
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN: regridded component vectors (3D arrays)
% OUTPUT:    
%   vel: regridded velocities - plate motion 
%   
% Andrew Watson     28-06-2022
%                                                                  
%=================================================================

plt = 0;

%% setup

% load plate motion vectors from file
plate_vels = readmatrix(par.plate_motion_file);

% crop to the overall region covered by the insar
insar_poly = [x(1) y(1); x(end) y(1); x(end) y(end); x(1) y(end)];
[in_poly,~] = inpolygon(plate_vels(:,1),plate_vels(:,2),insar_poly(:,1),insar_poly(:,2));
plate_vels_crop = plate_vels(in_poly,:);

% interpolate onto insar coords
[xx,yy] = meshgrid(x,y);
plate_E = griddata(plate_vels_crop(:,1),plate_vels_crop(:,2),plate_vels_crop(:,3),xx,yy);
plate_N = griddata(plate_vels_crop(:,1),plate_vels_crop(:,2),plate_vels_crop(:,4),xx,yy);

% for plotting
if plt == 1; vel_orig = vel; load('plotting/cpt/vik.mat'); end

%% loop through each velocity field

for ii = 1:size(vel,3)
    
    % project into line-of-sight
    plate_los = (plate_E .* compE(:,:,ii)) + (plate_N .* compN(:,:,ii));
    
    % set a common reference (find point closest to zero vel)
    [~,min_ind] = min(abs(vel(:,:,ii)),[],'all','linear','omitnan');
    [ref_yind,ref_xind] = ind2sub(size(vel(:,:,ii)),min_ind);    
    vel(:,:,ii) = vel(:,:,ii) - vel(ref_yind,ref_xind,ii);
    plate_los = plate_los - plate_los(ref_yind,ref_xind);
    
    % apply correction
    vel(:,:,ii) = vel(:,:,ii) - plate_los;
    
    % optional plotting
    if plt == 1

        % limits
        clim = [-10 10];
        [~,x_ind,y_ind] = crop_nans(vel(:,:,ii),x,y);
        lonlim = x([x_ind(1) x_ind(end)]); latlim = y([y_ind(1) y_ind(end)]);
        
        f = figure();
        f.Position([1 3 4]) = [600 1600 600];
        tiledlayout(1,3,'TileSpacing','compact')

        nexttile; hold on
        imagesc(x,y,vel_orig(:,:,ii),'AlphaData',~isnan(vel_orig(:,:,ii))); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(vik); caxis(clim)
        title('Original InSAR vel')

        nexttile; hold on
        imagesc(x,y,plate_los,'AlphaData',~isnan(plate_los)); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(vik); caxis(clim)
        title('Plate motion with shared ref')

        nexttile; hold on
        imagesc(x,y,vel(:,:,ii),'AlphaData',~isnan(vel(:,:,ii))); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(vik); caxis(clim)
        title('InSAR - plate')
        
    end
    
end





% [~,xmin_ind] = min(abs(plate_vels(:,1)-x(1)));
% [~,xmax_ind] = min(abs(plate_vels(:,1)-x(end)));
% [~,ymin_ind] = min(abs(plate_vels(:,2)-y(1)));
% [~,ymax_ind] = min(abs(plate_vels(:,2)-y(end)));
% plate_vels = plate_vels(ymin_ind:ymax_ind,xmin_ind:xmax_ind);