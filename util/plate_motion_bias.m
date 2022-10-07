function [vel] = plate_motion_bias(par,x,y,vel,compE,compN,asc_frames_ind,desc_frames_ind)
%=================================================================
% function plate_motion_bias()
%-----------------------------------------------------------------
% Mitigate the "reference frame effect" caused by rigid plate motions by
% subtracting No-Net-Rotation plate motion in ITRF from the LOS
% velocities.
%
% The input plate motions can be generated using the UNAVCO plate motion
% calculator, using ITRF2014 and No-Net-Rotation.
%
% The input plate motion vector file should contain (at least) the
% following columns in positions 1-4 : lon, lat, East vel, North vel.
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN: regridded component vectors (3D arrays)
%   asc_frames_ind, desc_frames_ind: indices in vel for asc and desc
%       frames/tracks
% OUTPUT:    
%   vel: regridded velocities - plate motion 
%   
% Andrew Watson     28-06-2022
%                                                                  
%=================================================================

%% setup

% load plate motion vectors from file
plate_vels = readmatrix(par.plate_motion_file);

% crop to the overall region covered by the insar, expanding the area by
% 20% of what the InSAR covers
buffer_x = (x(end)-x(1)).*0.2;
buffer_y = (y(end)-y(1)).*0.2;

insar_poly = [x(1)-buffer_x y(1)-buffer_y; 
                x(end)+buffer_x y(1)-buffer_y; 
                x(end)+buffer_x y(end)+buffer_y;
                x(1)-buffer_x y(end)+buffer_y];
[in_poly,~] = inpolygon(plate_vels(:,1),plate_vels(:,2),insar_poly(:,1),insar_poly(:,2));
plate_vels_crop = plate_vels(in_poly,:);

% interpolate onto insar coords
[xx,yy] = meshgrid(x,y);
plate_E = griddata(plate_vels_crop(:,1),plate_vels_crop(:,2),plate_vels_crop(:,3),xx,yy);
plate_N = griddata(plate_vels_crop(:,1),plate_vels_crop(:,2),plate_vels_crop(:,4),xx,yy);

% for plotting
if par.plt_plate_motion == 1; vel_orig = vel; load('plotting/cpt/vik.mat'); end

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
    if par.plt_plate_motion_indv == 1

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
    
    % report progress
    disp([num2str(ii) '/' num2str(size(vel,3)) ' complete'])
    
end

%% plot all corrected frames

if par.plt_plate_motion == 1
    
    % set plotting parameters
    lonlim = [min(x) max(x)];
    latlim = [min(y) max(y)];
    clim = [-10 10];
    load('plotting/cpt/vik.mat')
    
    % reload borders for ease
    if par.plt_borders == 1
        borders = load(par.borders_file);
    else
        borders = [];
    end
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'After plate motion correction')
    
    % plot ascending tracks
    t(1) = nexttile; hold on
    plt_data(x,y,vel(:,:,asc_frames_ind),lonlim,latlim,clim,'Ascending (mm/yr)',[],borders)
    colormap(t(1),vik)
    
    % plot descending tracks
    t(2) = nexttile; hold on
    plt_data(x,y,vel(:,:,desc_frames_ind),lonlim,latlim,clim,'Descending (mm/yr)',[],borders)
    colormap(t(2),vik)

end
