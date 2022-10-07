%% mask_multiple_frames.m
% Uses the .h5 output from licsbas to produce masks for licsbas velocities
% based on the standard licsbas noise indices.
% Written to be faster than rerunning LiCSBAS15_mask_ts.py for multiple
% frames, and allows for quick previews of coverage.
%
% Andrew Watson     19-07-2022

disp('Beginning run')

config_file = '/scratch/eearw/decomp_frame_vels/conf/iran_gacos.conf';

addpath ../plotting

%% toggles

plt_single_frame = 0;
plt_single_indices = 0;
new_masks = 1;
plt_new_masked_vels = 1;
plt_asc_desc_masks = 1;

%% setup

% coh_avg, n_unw, vstd, maxTlen, n_gap, stc, n_ifg_noloop, n_loop_err, resid_rms
thresholds = [0.2 0 20 3 0 100 200 75 3];

file_names = {'vel.geo.tif' 'coh_avg.geo.tif' 'n_unw.geo.tif' 'vstd.geo.tif' 'maxTlen.geo.tif' 'n_gap.geo.tif'...
     'stc.geo.tif' 'n_ifg_noloop.geo.tif' 'n_loop_err.geo.tif' 'resid_rms.geo.tif'};

%% read parameter file
% use the standard config file setup, assuming fixed file names

disp('Loading parameter file')

[par,insarpar] = readparfile(config_file);

%% load inputs

disp('Loading inputs')

% vel and noise indices stored as a 3d array in the cell array frame_data.
% order is : vel, coh_avg, n_unw, vstd, maxTlen, n_gap, stc, n_ifg_noloop, n_loop_err, resid_rms

% number of velocity maps inputted
nframes = length(insarpar.dir);

% pre-allocate
frames = cell(1,nframes);
lon = cell(1,nframes); lat = cell(size(lon));
% lon_comp = cell(size(lon)); lat_comp = cell(size(lon));
dx = cell(size(lon)); dy = cell(size(lon));

frame_data = cell(1,nframes);

% for each frame directory
for ii = 1:nframes
    
    disp(['Loading ' insarpar.dir{ii}])

    % extract the frame name
    frames(ii) = regexp(insarpar.dir{ii},'\d*[AD]_\d*_\d*','match');

    % read noise indices from h5 file
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vel '*']);
    [lon{ii},lat{ii},vel,dx{ii},dy{ii}] ...
        = read_geotiff([insarpar.dir{ii} namestruct.name]);

    % test that vel contains valid pixels, and remove if not
    if sum(~isnan(vel),'all') == 0
        disp([insarpar.dir{ii} ' is all nans - removing'])
        lon(ii) = []; lat(ii) = []; vel(ii) = []; dx(ii) = []; dy(ii) = [];
        continue
    end
    
    % write into 3d array with space for noise indices
    frame_data{ii} = zeros([size(vel) 10]);
    frame_data{ii}(:,:,1) = vel; 
    
    if new_masks == 1
        % load noise indices
        if isfile([insarpar.dir{ii} 'n_loop_err.geo.tif'])
            [~,~,frame_data{ii}(:,:,2),~,~] = read_geotiff([insarpar.dir{ii} 'coh_avg.geo.tif']);
            [~,~,frame_data{ii}(:,:,3),~,~] = read_geotiff([insarpar.dir{ii} 'n_unw.geo.tif']);
            [~,~,frame_data{ii}(:,:,4),~,~] = read_geotiff([insarpar.dir{ii} 'vstd.geo.tif']);
            [~,~,frame_data{ii}(:,:,5),~,~] = read_geotiff([insarpar.dir{ii} 'maxTlen.geo.tif']);
            [~,~,frame_data{ii}(:,:,6),~,~] = read_geotiff([insarpar.dir{ii} 'n_gap.geo.tif']);
            [~,~,frame_data{ii}(:,:,7),~,~] = read_geotiff([insarpar.dir{ii} 'stc.geo.tif']);
            [~,~,frame_data{ii}(:,:,8),~,~] = read_geotiff([insarpar.dir{ii} 'n_ifg_noloop.geo.tif']);
            [~,~,frame_data{ii}(:,:,9),~,~] = read_geotiff([insarpar.dir{ii} 'n_loop_err.geo.tif']);
            [~,~,frame_data{ii}(:,:,10),~,~] = read_geotiff([insarpar.dir{ii} 'resid_rms.geo.tif']);
        end

        % mask noise indices with unmasked vel (getting nans around frame)
        frame_data{ii}(isnan(repmat(vel,1,1,10))) = nan;
        
    elseif plt_single_indices ~= 0     
        if isfile([insarpar.dir{ii} 'n_loop_err.geo.tif'])
            [~,~,frame_data{ii}(:,:,plt_single_indices),~,~] = read_geotiff([insarpar.dir{ii} file_names{plt_single_indices}]);
            frame_data{ii}(isnan(repmat(vel,1,1,10))) = nan;
        end
        
    end
    
    % load mask
    if par.usemask == 1
        namestruct = dir([insarpar.dir{ii} '*' insarpar.id_mask '*']);
        [~,~,mask{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name]);
    end
    
end

% get indices of ascending and descending frames
asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), frames));
desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), frames));

% fault traces
if par.plt_faults == 1
    fault_trace = readmatrix(par.faults_file);
else
    fault_trace = [];
end

% gnss vels
load(par.gnss_file);

% borders
if par.plt_borders == 1
    borders = load(par.borders_file);
else
    borders = [];
end

% colour palettes
load('../plotting/cpt/vik.mat')
load('../plotting/cpt/batlow.mat')

%% generate new masks

if new_masks == 1
    for ii = 1:nframes

        % reset mask
        mask{ii} = ones(size(frame_data{ii}(:,:,1)));
        mask{ii}(isnan(frame_data{ii}(:,:,1))) = nan;

        % apply new thresholds
        mask{ii}(frame_data{ii}(:,:,2) < thresholds(1)) = 0;
        mask{ii}(frame_data{ii}(:,:,3) < thresholds(2)) = 0;
        mask{ii}(frame_data{ii}(:,:,4) > thresholds(3)) = 0;
        mask{ii}(frame_data{ii}(:,:,5) < thresholds(4)) = 0;
        mask{ii}(frame_data{ii}(:,:,6) > thresholds(5)) = 0;
        mask{ii}(frame_data{ii}(:,:,7) > thresholds(6)) = 0;
        mask{ii}(frame_data{ii}(:,:,8) > thresholds(7)) = 0;
        mask{ii}(frame_data{ii}(:,:,9) > thresholds(8)) = 0;
        mask{ii}(frame_data{ii}(:,:,10) > thresholds(9)) = 0;

    end
end

%% plot single frame

if plt_single_frame ~= 0
    plot_frame_indices(frame_data{plt_single_frame},mask{plt_single_frame},thresholds,vik)
end

%% plot single indices

if plt_single_indices ~= 0
     
    clim = [0 100];
    
    % asc
    [p,~] = numSubplots(length(asc_frames_ind));
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    tiledlayout(p(1),p(2),'TileSpacing','compact')
    
    for ii = 1:length(asc_frames_ind)       
        nexttile
        plt_data(lon{asc_frames_ind(ii)},lat{asc_frames_ind(ii)},...
            frame_data{asc_frames_ind(ii)}(:,:,plt_single_indices),[],[],clim,frames{asc_frames_ind(ii)},[],[])
%         colormap(t(1),vik)
    end
    
    
    % desc
    [p,~] = numSubplots(length(desc_frames_ind));
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    tiledlayout(p(1),p(2),'TileSpacing','compact')
    
    for ii = 1:length(desc_frames_ind)       
        nexttile
        plt_data(lon{desc_frames_ind(ii)},lat{desc_frames_ind(ii)},...
            frame_data{desc_frames_ind(ii)}(:,:,plt_single_indices),[],[],clim,frames{desc_frames_ind(ii)},[],[])
%         colormap(t(1),vik)
    end
    
    
end

%% plot multiple

if plt_new_masked_vels == 1

    % set plotting parameters
    lonlim = [min(cellfun(@min,lon)) max(cellfun(@max,lon))];
    latlim = [min(cellfun(@min,lat)) max(cellfun(@max,lat))];
    clim = [-10 10];

    % temp apply mask
    vel_tmp = cell(1,nframes);
    for ii = 1:nframes
        vel_tmp{ii} = frame_data{ii}(:,:,1);
        vel_tmp{ii}(mask{ii}==0) = nan;
        vel_tmp{ii} = vel_tmp{ii} - median(vel_tmp{ii},'all','omitnan');
    end

    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    tiledlayout(1,2,'TileSpacing','compact')

    % plot ascending tracks
    t(1) = nexttile; hold on
    plt_data(lon(asc_frames_ind),lat(asc_frames_ind),vel_tmp(asc_frames_ind),...
        lonlim,latlim,clim,'Ascending (mm/yr)',[],borders)
    colormap(t(1),vik)

    % plot descending tracks
    t(2) = nexttile; hold on
    plt_data(lon(desc_frames_ind),lat(desc_frames_ind),vel_tmp(desc_frames_ind),...
        lonlim,latlim,clim,'Descending (mm/yr)',[],borders)
    colormap(t(2),vik)

    clear vel_tmp

end

%% functions ==============================================================

function plot_frame_indices(data,mask,thresholds,vik)
% 3D array order: vel, coh_avg, n_unw, vstd, maxTlen, n_gap, stc, n_ifg_noloop, n_loop_err, resid_rms

figure()
tiledlayout(3,4,'TileSpacing','compact')

nexttile
imagesc(data(:,:,1),'AlphaData',~isnan(data(:,:,1)))
title('vel')
colorbar
caxis([-20 20])
ax = gca; colormap(ax,vik)

nexttile
masked_vel = data(:,:,1);
masked_vel(mask==0) = nan;
imagesc(masked_vel,'AlphaData',~isnan(masked_vel))
title('vel\_mask')
colorbar
caxis([-20 20])
ax = gca; colormap(ax,vik)

nexttile
imagesc(mask,'AlphaData',~isnan(mask))
title('mask')
colorbar

nexttile
imagesc(data(:,:,2),'AlphaData',~isnan(data(:,:,2)))
title(['coh\_avg (' num2str(thresholds(1)) ')'])
colorbar

nexttile
imagesc(data(:,:,3),'AlphaData',~isnan(data(:,:,3)))
title(['n\_unw (' num2str(thresholds(2)) ')'])
colorbar

nexttile
imagesc(data(:,:,4),'AlphaData',~isnan(data(:,:,4)))
title(['vstd (' num2str(thresholds(3)) ')'])
colorbar

nexttile
imagesc(data(:,:,5),'AlphaData',~isnan(data(:,:,5)))
title(['maxTlen (' num2str(thresholds(4)) ')'])
colorbar

nexttile
imagesc(data(:,:,6),'AlphaData',~isnan(data(:,:,6)))
title(['n\_gap (' num2str(thresholds(5)) ')'])
colorbar

nexttile
imagesc(data(:,:,7),'AlphaData',~isnan(data(:,:,7)))
title(['stc (' num2str(thresholds(6)) ')'])
colorbar

nexttile
imagesc(data(:,:,8),'AlphaData',~isnan(data(:,:,8)))
title(['n\_ifg\_noloop (' num2str(thresholds(7)) ')'])
colorbar

nexttile
imagesc(data(:,:,9),'AlphaData',~isnan(data(:,:,9)))
title(['n\_loop\_err (' num2str(thresholds(8)) ')'])
colorbar

nexttile
imagesc(data(:,:,10),'AlphaData',~isnan(data(:,:,10)))
title(['resid\_rms (' num2str(thresholds(9)) ')'])
colorbar

end