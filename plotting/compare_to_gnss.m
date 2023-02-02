%% compare_to_gnss.m
% Plot the residual between gnss velocities and decomposed or line-of-sight
% insar velocities.
% GNSS files is expected to have lon, lat, east, north in first 4 cols.
% Compares to the median of the pixels within dist_threshold of each gnss
% station.
%
% Andrew Watson     2022-09-12

addpath ../util/

%% setup

% direction of velocities ['east' 'north' 'los']
vel_direction = 'north';

% distance threhold, in same coords as velocities (likely degrees)
dist_or_nearest = 0; %0=dist, 1=nearest
dist_threshold = 0.10;

% fit some function to the gnss residuals
fit_func_to_resid = 0;
remove_offset = 0;

%% load

% track = '006D';
% in_dir = '/scratch/eearw/iran_frame_vels/long_track_comparisons/plate_corrd_2km/';

% component files (required for los)
% compE_file = [in_dir track '_compE.geo.tif'];
% compN_file = [in_dir track '_compN.geo.tif'];

% vels
% vel_file = [in_dir track '_vel.geo.tif'];
% vel_file = '/scratch/eearw/decomp_frame_vels/out/2km_for_plotting/iran_gacos_ml2_vE.geo.tif';
% vel_file = '/scratch/eearw/decomp_frame_vels/out/interp_test/iran_gacos_2km_vE.geo.tif';
% vel_file = '/scratch/eearw/kriging/out/khor_cleaned_extended/testing/vn_interpolated_synth2.mat';
vel_file = '/nfs/a285/homes/eearw/gmt/thesis/chp4/gnss_kriging/data/mkdf_pinning/vn_interpolated_mkdf1.mat';

% change load method for tif and mat files
[~,~,ext] = fileparts(vel_file);
if ext == '.tif'
    [lon,lat,vel,~,~] = read_geotiff(vel_file);
elseif ext == '.mat'
    load(vel_file)
end

% gnss
% gnss_file = '/scratch/eearw/decomp_frame_vels/gnss/khor/cleaned_stations/khor_vert_10mm_gf7_buff01.csv';
gnss_file = '/scratch/eearw/decomp_frame_vels/gnss/khor/cleaned_stations/khor_vert_10mm_gf7_buff01_withMkdf1.csv';
gnss = readmatrix(gnss_file);

% borders for plotting
borders = load('/nfs/a285/homes/eearw/velmap/plotting/borderdata.mat');

%% tidying

% crop padding nans on vel
[vel,crop_xind,crop_yind,lon,lat] = crop_nans(vel,lon,lat);

% remove any gnss vels not within the area of the vel (including nans)
outside_area = (gnss(:,1) < min(lon) | gnss(:,1) > max(lon)) | ...
    (gnss(:,2) < min(lat) | gnss(:,2) > max(lat));
gnss(outside_area,:) = [];

%% select gnss component

switch vel_direction
    case 'east'
        gnss_vel = gnss(:,[1 2 3]);
        
    case 'north'
        gnss_vel = gnss(:,[1 2 4]);
        
    case 'los'
        % pre-al
        gnss_vel = nan(size(gnss,1),3);
        gnss_vel(:,1:2) = gnss(:,1:2);
        
        % load component files
        [~,~,compE,~,~] = read_geotiff(compE_file);
        [~,~,compN,~,~] = read_geotiff(compN_file);
        
        % crop
        compE = compE(crop_yind,crop_xind);
        compN = compN(crop_yind,crop_xind);
        
        % for each gnss, find nearest point in compE and compN, and project
        % into los
        for ii = 1:size(gnss,1)

            [~,ind_x] = min(abs(lon-gnss_vel(ii,1)));
            [~,ind_y] = min(abs(lat-gnss_vel(ii,2)));
            
            gnss_vel(ii,3) = gnss(ii,3).*compE(ind_y,ind_x) ...
                + gnss(ii,4).*compN(ind_y,ind_x);
            
        end
        
end

%% calc resid

% pre-al
resid = nan(size(gnss_vel,1),3);
resid(:,1:2) = gnss_vel(:,1:2);

% coords grid
[xx,yy] = meshgrid(lon,lat);

% loop through gnss
for ii = 1:size(gnss,1)
      
    if dist_or_nearest == 0
        % distance from gnss
        dist_from_gnss = sqrt((xx-gnss_vel(ii,1)).^2 + (yy-gnss_vel(ii,2)).^2);

        % residual
        resid(ii,3) = gnss_vel(ii,3) - median(vel(dist_from_gnss<=dist_threshold),'omitnan');
        
        
    elseif dist_or_nearest == 1
        % nearest point
        [~,ind_x] = min(abs(lon-gnss_vel(ii,1)));
        [~,ind_y] = min(abs(lat-gnss_vel(ii,2)));
        
        resid(ii,3) = gnss_vel(ii,3) - vel(ind_y,ind_x);
        
    end
    
end

% clear nans (where gnss and vel don't overlap)
resid(isnan(resid(:,3)),:) = [];

%% remove offset

if remove_offset == 1
    
%     vel = vel + mean(resid(:,3));
    gnss_vel(:,3) = gnss_vel(:,3) - mean(resid(:,3));
    
    % loop through gnss
    for ii = 1:size(gnss,1)

        if dist_or_nearest == 0
            % distance from gnss
            dist_from_gnss = sqrt((xx-gnss_vel(ii,1)).^2 + (yy-gnss_vel(ii,2)).^2);

            % residual
            resid(ii,3) = gnss_vel(ii,3) - median(vel(dist_from_gnss<=dist_threshold),'omitnan');


        elseif dist_or_nearest == 1
            % nearest point
            [~,ind_x] = min(abs(lon-gnss_vel(ii,1)));
            [~,ind_y] = min(abs(lat-gnss_vel(ii,2)));

            resid(ii,3) = gnss_vel(ii,3) - vel(ind_y,ind_x);

        end

    end
      
end

% clear nans (where gnss and vel don't overlap)
resid(isnan(resid(:,3)),:) = [];

%% fit to residuals

if fit_func_to_resid == 1
    
    % fit polynomial surface
    fitobject = fit([resid(:,1) resid(:,2)],resid(:,3),'poly11');

    % evaluate for all points
    [xx,yy] = meshgrid(lon,lat);
    resid_surface = fitobject(xx,yy);
    
    % residual with stations
    resid_surface_gnss = resid(:,3) - fitobject(resid(:,1),resid(:,2));
    
    % plot
    figure()
    tiledlayout(2,1,'TileSpacing','compact')
    
    nexttile(); hold on
    imagesc(lon,lat,resid_surface)
    scatter(resid(:,1),resid(:,2),40,resid_surface_gnss,'Filled','MarkerEdgeColor','k')
    axis xy
    colorbar
    title('Fit to GNSS residuals')
    
    nexttile(); hold on
    histogram(resid_surface_gnss,20);
    
    
    load('cpt/vik.mat')
    figure()
    imagesc(lon,lat,vel+resid_surface,'AlphaData',~isnan(vel))
    colorbar; colormap(vik); caxis([-10 10])
    xlim(lonlim); ylim(latlim)
    
end

%% plot

load('cpt/vik.mat')

lonlim = [min(lon) max(lon)];
latlim = [min(lat) max(lat)];

% figure(); hold on
% tiledlayout(2,1,'TileSpacing','compact')
% 
% nexttile(); hold on
% imagesc(lon,lat,vel,'AlphaData',~isnan(vel))
% scatter(resid(:,1),resid(:,2),40,resid(:,3),'Filled','MarkerEdgeColor','k')
% colorbar; colormap(vik); caxis([-5 5])
% xlim(lonlim); ylim(latlim)
% 
% nexttile()
% histogram(resid(:,3),20)
% ylabel('Residual (mm/yr)')

f = figure(); hold on
% f.Position = [50 700 600 1000];
imagesc(lon,lat,vel,'AlphaData',~isnan(vel))
scatter(resid(:,1),resid(:,2),70,resid(:,3),'Filled','MarkerEdgeColor','k')
for ii = 1:length(borders.places); plot(borders.lon{ii},borders.lat{ii},'k'); end
colorbar; colormap(vik); caxis([-5 5])
xlim(lonlim); ylim(latlim)
% title(track)

f = figure();
% f.Position = [50 100 600 500];
histogram(resid(:,3),20)
title(['mean = ' num2str(mean(resid(:,3))) ', SD = ' num2str(std(resid(:,3)))])
xlabel('Residual (mm/yr)')
