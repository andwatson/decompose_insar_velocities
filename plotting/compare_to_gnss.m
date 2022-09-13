%% compare_to_gnss.m
% Plot residual between gnss velocities and decomposed insar velocities.
% Currently hardcoded to compare to East GNSS vel.
% Compares to the median of the pixels within dist_threshold of each gnss
% station.
%
% Andrew Watson     2022-09-12

dist_threshold = 0.1;

%% load

% vels
vel_file = '/scratch/eearw/decomp_frame_vels/out/2km_for_plotting/iran_gacos_ml2_vE.geo.tif';
[lon,lat,vel,~,~] = read_geotiff(vel_file);

% gnss
gnss_file = '/scratch/eearw/decomp_frame_vels/gnss/khor/cleaned_stations/khor_vert_10mm_gf7_buff01.csv';
gnss = readmatrix(gnss_file);

%% calc resid

% pre-al
resid = nan(size(gnss,1),3);
resid(:,1:2) = gnss(:,1:2);

% coords grid
[xx,yy] = meshgrid(lon,lat);

% loop through gnss
for ii = 1:size(gnss,1)
    
    % distance from gnss
    dist_from_gnss = sqrt((xx-gnss(ii,1)).^2 + (yy-gnss(ii,2)).^2);
    
    % residual
    resid(ii,3) = gnss(ii,3) - median(vel(dist_from_gnss<=dist_threshold),'omitnan');
    
end

% clear nans (where gnss and vel don't overlap)
resid(isnan(resid(:,3)),:) = [];

%% plot

load('cpt/vik.mat')

lonlim = [min(lon) max(lon)];
latlim = [min(lat) max(lat)];

figure(); hold on
tiledlayout(2,1,'TileSpacing','compact')

nexttile(); hold on
imagesc(lon,lat,vel,'AlphaData',~isnan(vel))
scatter(resid(:,1),resid(:,2),40,resid(:,3),'Filled','MarkerEdgeColor','k')
colorbar; colormap(vik); caxis([-5 5])
xlim(lonlim); ylim(latlim)

nexttile()
histogram(resid(:,3))
ylabel('Residual (mm/yr)')