%% compare_decomp_tifs.m
% Compare two velocity fields stored as tifs, assumed to be the same
% component and on the same grid.
%
% Andrew Watson     2022-09-13

%% load inputs

vel_file1 = '/scratch/eearw/decomp_frame_vels/out/thesis/20230202/iran_gacos_ml1_vE.geo.tif';
vel_file2 = '/scratch/eearw/decomp_frame_vels/out/thesis/vEvU/iran_gacos_ml1_vE.geo.tif';

[lon,lat,vel1,~,~] = read_geotiff(vel_file1);
[~,~,vel2,~,~] = read_geotiff(vel_file2);

% for plotting
load('cpt/vik.mat')
borders = load('/nfs/a285/homes/eearw/velmap/plotting/borderdata.mat');

%% calculate residual and plot

residual = vel1 - vel2;
resid_vec = residual(~isnan(residual));

lonlim = [min(lon) max(lon)];
latlim = [min(lat) max(lat)];
clim = [-10 10];

f = figure();
f.Position([1 3 4]) = [600 1600 1200];
tiledlayout(1,3,'TileSpacing','compact')

t(1) = nexttile; hold on
plt_data(lon,lat,vel1,lonlim,latlim,clim,'vel1 (mm/yr)',[],borders)
colormap(t(1),vik)

t(2) = nexttile; hold on
plt_data(lon,lat,vel2,lonlim,latlim,clim,'vel2 (mm/yr)',[],borders)
colormap(t(2),vik)

t(3) = nexttile; hold on
plt_data(lon,lat,residual,lonlim,latlim,clim,'vel1 - vel2 (mm/yr)',[],borders)
colormap(t(3),vik)