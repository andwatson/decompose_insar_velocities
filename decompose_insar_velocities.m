%% decompose_insar_velocities.m
%
% A continuation of the original merge_all_NEU scripts, designed to work
% with any number of input frames (not just the original 5 mrf ones).
%
% The setup is now done with a parameter file for flexibility.
%
% Decomposes line-of-sight velocities from LiCSBAS into a combination of
% North, East, and Up velocities, assuming no covariance between adjacent
% pixels. All frames are interpolated onto a new grid.
%
% InSAR velocities may be referenced to a GNSS velocity field (tie2gnss).
% The North GNSS velocities are subtracted from the InSAR velocities before
% the decomposition is performed.
%
% GNSS file should be a .mat file containing structure with the following:
%   gnss_field.x - vector of x axis coords (n)
%   gnss_field.y - vector of y axis coords (m)
%   gnss_field.N - grid of north gnss vels (mxn)
%   gnss_field.E - grid of east gnss vels (mxn)
% Optional extras:
%   gnss_field.sN - grid of north 1-sigma uncertainties (mxn)
%   gnss_field.sE - grid of east 1-sigma uncertainties (mxn)
%
% Andrew Watson     26-04-2021

disp('Beginning run')

config_file = '/scratch/eearw/decomp_frame_vels/conf/thesis_ds.conf';

% add subdirectory paths
addpath util plotting

% begin timer
tic

%% read parameter file

disp('Loading parameter file')

[par,insarpar] = readparfile(config_file);

%% check that inputs exist

disp('Checking that data exists')

to_remove = false(1,length(insarpar.dir));

for ii = 1:length(insarpar.dir)
    
    % get name of velocity file
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vel '*']);
    
    % test existance, record if config file if missing
    if isempty(namestruct)   
        disp(['Removing ' insarpar.dir{ii}])
        to_remove(ii) = true;        
    end
    
end

% remove missing file dirs
insarpar.dir(to_remove) = [];

%% load inputs

disp('Loading inputs')

% number of velocity maps inputted
nframes = length(insarpar.dir);

% pre-allocate
frames = cell(1,nframes);
lon = cell(1,nframes); lat = cell(size(lon));
lon_comp = cell(size(lon)); lat_comp = cell(size(lon));
dx = cell(size(lon)); dy = cell(size(lon));
vel = cell(size(lon)); vstd = cell(size(lon)); mask = cell(size(lon));
compE = cell(size(lon)); compN = cell(size(lon)); compU = cell(size(lon));

% for each velocity map
for ii = 1:nframes
    
    disp(['Loading ' insarpar.dir{ii}])

    % extract the frame name
    frames(ii) = regexp(insarpar.dir{ii},'\d*[AD]_\d*_\d*','match');

    % load velocities
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vel '*']);
    [lon{ii},lat{ii},vel{ii},dx{ii},dy{ii}] ...
        = read_geotiff([insarpar.dir{ii} namestruct.name]);

    % test that vel contains valid pixels, and remove if not
    if sum(~isnan(vel{ii}),'all') == 0
        disp([insarpar.dir{ii} ' is all nans - removing'])
        lon(ii) = []; lat(ii) = []; vel(ii) = []; dx(ii) = []; dy(ii) = [];
        continue
    end

    % load velocity errors
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vstd '*']);
    [~,~,vstd{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name]);

    % load East, North, and Up components
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_e '*']);
    [lon_comp{ii},lat_comp{ii},compE{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name]);

    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_n '*']);
    [~,~,compN{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name]);

    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_u '*']);
    [~,~,compU{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name]);
    
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
load('plotting/cpt/vik.mat')
load('plotting/cpt/batlow.mat')

%% preview inputs

if par.plt_input_vels == 1
    
    % set plotting parameters
    lonlim = [min(cellfun(@min,lon)) max(cellfun(@max,lon))]; 
    latlim = [min(cellfun(@min,lat)) max(cellfun(@max,lat))]; 
    clim = [-10 10];
    
    % temp apply mask
    vel_tmp = vel;
    for ii = 1:nframes
        vel_tmp{ii}(mask{ii}==0) = nan;
        vel_tmp{ii} = vel_tmp{ii} - median(vel_tmp{ii},'all','omitnan');
    end
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'Input velocities')
    
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

%% downsample unit vectors if required
% ENU from licsbas may not be downsampled to the same level as vel by
% default, so check the sizes and downsample if neccessary. 
disp('Checking if unit vectors need downsampling')

for ii = 1:nframes
    
    % calculate downsample factor from matrix size
    comp_dsfac = round(size(compE{ii},1)./size(vel{ii},1));
    
    if comp_dsfac > 1
        
        disp(['Downsampling unit vectors : ' insarpar.dir{ii}])
    
        % downsample components
        [compE{ii},lon_comp{ii},lat_comp{ii}] = downsample_array(compE{ii},...
            comp_dsfac,comp_dsfac,'mean',lon_comp{ii},lat_comp{ii});
        [compN{ii},~,~] = downsample_array(compN{ii},comp_dsfac,comp_dsfac,'mean');
        [compU{ii},~,~] = downsample_array(compU{ii},comp_dsfac,comp_dsfac,'mean');
        
    end

end

%% downsample

if par.ds_factor > 0
    
    disp(['Downsampling by a factor of ' num2str(par.ds_factor)])
    
    for ii = 1:nframes
    
        [vel{ii},lon{ii},lat{ii}] ...
            = downsample_array(vel{ii},par.ds_factor,par.ds_factor,par.ds_method,lon{ii},lat{ii});        
        [vstd{ii},~,~] = downsample_array(vstd{ii},par.ds_factor,par.ds_factor,par.ds_method);
        
        [compE{ii},lon_comp{ii},lat_comp{ii}] ...
            = downsample_array(compE{ii},par.ds_factor,par.ds_factor,par.ds_method,lon_comp{ii},lat_comp{ii});
        [compN{ii},~,~] = downsample_array(compN{ii},par.ds_factor,par.ds_factor,par.ds_method);
        [compU{ii},~,~] = downsample_array(compU{ii},par.ds_factor,par.ds_factor,par.ds_method);
        
        if par.usemask == 1
            [mask{ii},~,~] ...
                =  downsample_array(mask{ii},par.ds_factor,par.ds_factor,par.ds_method);
        end
        
        if (mod(ii,round(nframes./10))) == 0
            disp([num2str(round((ii./nframes)*100)) '% completed']);
        end
        
    end
    
end

%% unify grids

disp('Unifying grids')

% limits and intervals for new grid
x_regrid = min(cellfun(@min,lon)) : min(cellfun(@min,dx)) : max(cellfun(@max,lon));
y_regrid = min(cellfun(@min,lat)) : min(cellfun(@min,dy)) : max(cellfun(@max,lat));
dx = min(cellfun(@min,dx)); dy = min(cellfun(@min,dy));

% new grid to project data onto
[xx_regrid,yy_regrid] = meshgrid(x_regrid,y_regrid); 

% pre-allocate
vel_regrid = zeros([size(xx_regrid) nframes]);
vstd_regrid = zeros(size(vel_regrid)); mask_regrid = zeros(size(vel_regrid));
compE_regrid = zeros(size(vel_regrid)); compN_regrid = zeros(size(vel_regrid));
compU_regrid = zeros(size(vel_regrid));

for ii = 1:nframes
    
    % reduce interpolation area using coords of frame
    [~,xind_min] = min(abs(x_regrid-lon{ii}(1)));
    [~,xind_max] = min(abs(x_regrid-lon{ii}(end)));
    [~,yind_min] = min(abs(y_regrid-lat{ii}(end)));
    [~,yind_max] = min(abs(y_regrid-lat{ii}(1)));    
    
    % interpolate vel, vstd, and optional mask
    [xx,yy] = meshgrid(lon{ii},lat{ii});
    vel_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
        = interp2(xx,yy,vel{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
        yy_regrid(yind_min:yind_max,xind_min:xind_max));
    
    vstd_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
        = interp2(xx,yy,vstd{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
        yy_regrid(yind_min:yind_max,xind_min:xind_max));
    
    if par.usemask == 1
        mask_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
            = interp2(xx,yy,mask{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
            yy_regrid(yind_min:yind_max,xind_min:xind_max));
    end
    
    % ENU may not be downsampled by licsbas, hence different xx yy grids
    [xx,yy] = meshgrid(lon_comp{ii},lat_comp{ii});
    compE_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
        = interp2(xx,yy,compE{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
        yy_regrid(yind_min:yind_max,xind_min:xind_max));
    compN_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
        = interp2(xx,yy,compN{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
        yy_regrid(yind_min:yind_max,xind_min:xind_max));
    compU_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
        = interp2(xx,yy,compU{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
        yy_regrid(yind_min:yind_max,xind_min:xind_max));
      
    if (mod(ii,round(nframes./10))) == 0
        disp([num2str(round((ii./nframes)*100)) '% completed']);
    end
    
end

if par.tie2gnss ~= 0
    [xx_gnss,yy_gnss] = meshgrid(gnss_field.x,gnss_field.y);
    gnss_E = interp2(xx_gnss,yy_gnss,gnss_field.E,xx_regrid,yy_regrid);
    gnss_N = interp2(xx_gnss,yy_gnss,gnss_field.N,xx_regrid,yy_regrid);
    gnss_sE = interp2(xx_gnss,yy_gnss,gnss_field.sE,xx_regrid,yy_regrid);
    gnss_sN = interp2(xx_gnss,yy_gnss,gnss_field.sN,xx_regrid,yy_regrid);
end

% change zeros to nans
vstd_regrid(vstd_regrid==0) = nan;
compE_regrid(compE_regrid==0) = nan;
compN_regrid(compN_regrid==0) = nan;
compU_regrid(compU_regrid==0) = nan;

% clear ungridded data
clear vel vstd compE compN compU

%% apply mask

if par.usemask == 1
    
    disp('Applying mask')
    
    mask_regrid(isnan(mask_regrid)) = 0;
    
    % account for previous downsampling
    mask_regrid(mask_regrid>=0.5) = 1; mask_regrid(mask_regrid<0.5) = 0;
    
    vel_regrid(mask_regrid==0) = NaN;
    
%     for ii=1:nframes
%         vel{ii}(mask{ii}==0) = nan;
%     end

end

%% plot ascending and descending masks
% Produces two masks that are useful for plotting.
% For overlaps, a given point is considered unmasked if it is unmasked in
% at least one of the masks.

if par.plt_mask_asc_desc == 1
    
    % look direction indices
    asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), frames));
    desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), frames));
    
    % sum masks
    mask_asc = sum(mask_regrid(:,:,asc_frames_ind),3);
    mask_desc = sum(mask_regrid(:,:,desc_frames_ind),3);
    
    % return to ones and zeros
    mask_asc(mask_asc>=1) = 1;
    mask_desc(mask_desc>=1) = 1;
    mask_asc(mask_asc<1) = 0;
    mask_desc(mask_desc<1) = 0;
    mask_asc(isnan(mask_asc)) = 0;
    mask_desc(isnan(mask_desc)) = 0;
    
    % plot
    lonlim = [min(x_regrid) max(x_regrid)];
    latlim = [min(y_regrid) max(y_regrid)];
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title('Combined masks for ascending and descending')
    
    t(1) = nexttile; hold on
    plt_data(x_regrid,y_regrid,mask_asc,lonlim,latlim,[],'Ascending mask',fault_trace,borders)
    
    t(2) = nexttile; hold on
    plt_data(x_regrid,y_regrid,mask_desc,lonlim,latlim,[],'Descending mask',fault_trace,borders)
    
end

%% merge frames along-track (and across-track)

if par.merge_tracks_along > 0
    
    disp('Merging frames along-track')
    
    [vel_regrid,compE_regrid,compN_regrid,compU_regrid,vstd_regrid,tracks] ...
        = merge_frames_along_track(par,x_regrid,y_regrid,vel_regrid,...
        frames,compE_regrid,compN_regrid,compU_regrid,vstd_regrid);
    
    % update number of frames and indexes if frames have been merged
    if par.merge_tracks_along == 2
        nframes = size(vel_regrid,3);
        asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), tracks));
        desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), tracks));
    end
    
    % run across-track merging
    if par.merge_tracks_across > 0        
        disp('Merging frame across-track. Not used in decomposition.')       
        merge_frames_across_track(par,x_regrid,y_regrid,vel_regrid,tracks,...
            compE_regrid,compU_regrid,vstd_regrid)
    end
    
end

%% reference frame bias correction
% remove "reference frame bias" caused by rigid plate motions.

if par.plate_motion == 1    
    [vel_regrid] = plate_motion_bias(par,x_regrid,y_regrid,vel_regrid,...
        compE_regrid,compN_regrid,asc_frames_ind,desc_frames_ind);
end

%% tie to gnss
% Shift InSAR velocities into the same reference frame as the GNSS
% velocities. Method is given by par.tie2gnss. 

if par.tie2gnss ~= 0
    [vel_regrid] = ref_to_gnss(par,xx_regrid,yy_regrid,...
        vel_regrid,compE_regrid,compN_regrid,gnss_E,gnss_N,frames);    
end

% calculate frame overlaps if requested
if par.frame_overlaps == 1
    disp('Calculating frame overlap statistics')
    frame_overlap_stats(vel_regrid,frames,compU_regrid);
end

%% identify pixels without both asc and desc coverage

asc_coverage = any(~isnan(vel_regrid(:,:,asc_frames_ind)),3);
desc_coverage = any(~isnan(vel_regrid(:,:,desc_frames_ind)),3);
both_coverage = all(cat(3,asc_coverage,desc_coverage),3);

%% velocity decomposition

disp('Inverting for E and U')

% check that at least one point has at least two look directions
if sum(both_coverage(:)) == 0
    error('No pixels with multiple look directions - did you provide more than one track?')
end

if par.decomp_method == 0 || par.decomp_method == 1
    
    % either remove gnss N, and decompose into vE and vU, or include gnss N
    % in the linear problem, and solve for vE, vU, and vN
    [m_east,m_up,var_east,var_up,condG_threshold_mask,var_threshold_mask] ...
        = vel_decomp(par,vel_regrid,vstd_regrid,compE_regrid,compN_regrid,...
        compU_regrid,gnss_N,gnss_sN,both_coverage);

elseif par.decomp_method == 2
    
    % decompose into vE and vUN, then split vUN into vU and vN (Qi's
    % method)
    [m_east,m_up,var_east,var_up,condG_threshold_mask,var_threshold_mask] ...
        = vel_decomp_vE_vUN(par,vel_regrid,vstd_regrid,compE_regrid,compN_regrid,...
        compU_regrid,gnss_N,gnss_sN,both_coverage);    
end
    
% [m_east,m_up,var_east,var_up,condG_threshold_mask,var_threshold_mask] ...
%     = null_line_decomp(par,vel_regrid,vstd_regrid,compE_regrid,compN_regrid,compU_regrid,both_coverage,asc_frames_ind,desc_frames_ind);

%% plot output velocities

lonlim = [min(x_regrid) max(x_regrid)];
latlim = [min(y_regrid) max(y_regrid)];
clim = [-10 10];

% East and vertical
f = figure();
f.Position([1 3 4]) = [600 2800 1200];
t = tiledlayout(1,2,'TileSpacing','compact');
title(t,'Decomposed velocities')

t(1) = nexttile; hold on
plt_data(x_regrid,y_regrid,m_up,lonlim,latlim,clim,'Vertical (mm/yr)',fault_trace,borders)
colormap(t(1),vik)

t(2) = nexttile; hold on
plt_data(x_regrid,y_regrid,m_east,lonlim,latlim,clim,'East (mm/yr)',fault_trace,borders)
colormap(t(2),vik)

% North and  coverage
f = figure();
f.Position([1 3 4]) = [600 2800 1200];
t = tiledlayout(1,2,'TileSpacing','compact');
title(t,'Decomposed velocities and coverage')

t(3) = nexttile; hold on
plt_data(x_regrid,y_regrid,gnss_N,lonlim,latlim,[],'North (mm/yr)',fault_trace,borders)

t(4) = nexttile; hold on
coverage = sum(~isnan(vel_regrid),3);
plt_data(x_regrid,y_regrid,coverage,lonlim,latlim,[],'Data coverage',fault_trace,borders)

%% plot velocity uncertainties

if par.plt_decomp_uncer == 1

    lonlim = [min(x_regrid) max(x_regrid)];
    latlim = [min(y_regrid) max(y_regrid)];
    clim = [0 4];

    f = figure();
    f.Position([1 3 4]) = [600 1600 800];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'Decomposed velocity uncertainties')

    t(1) = nexttile; hold on
    plt_data(x_regrid,y_regrid,var_up,lonlim,latlim,clim,'Vertical (mm/yr)',fault_trace,borders)
    colormap(t(1),batlow)

    t(2) = nexttile; hold on
    plt_data(x_regrid,y_regrid,var_east,lonlim,latlim,clim,'East (mm/yr)',fault_trace,borders)
    colormap(t(2),batlow)

end

%% plot variance and cond(G) threshold masks if used

if par.plt_threshold_masks == 1
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'Threshold masks')
    
    t(1) = nexttile; hold on
    plt_data(x_regrid,y_regrid,condG_threshold_mask,lonlim,latlim,[],'cond(G) mask',fault_trace,borders)
    
    t(2) = nexttile; hold on
    plt_data(x_regrid,y_regrid,var_threshold_mask,lonlim,latlim,[],'variance mask',fault_trace,borders)
    
end

%% save outputs

if par.save_geotif == 1
    
    % create georeference
    georef = georefpostings([min(y_regrid) max(y_regrid)],...
        [min(x_regrid) max(x_regrid)],size(m_up),'ColumnsStartFrom','south',...
        'RowsStartFrom','west');
    
    % write geotifs
    geotiffwrite([par.out_path par.out_prefix '_vU.geo.tif'],m_up,georef)
    geotiffwrite([par.out_path par.out_prefix '_vE.geo.tif'],m_east,georef)
    geotiffwrite([par.out_path par.out_prefix '_vN.geo.tif'],m_north,georef)
    
end

%% end

disp('Run complete')
toc