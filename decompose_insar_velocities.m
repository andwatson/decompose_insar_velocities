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
%
% Andrew Watson     26-04-2021

disp('Beginning run')

config_file = '/scratch/eearw/decomp_frame_vels/conf/iran_gacos.conf';

% add subdirectory paths
addpath util

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

%% unify grids

disp('Unifying grids')

% limits and intervals for new grid
x_regrid = min(cellfun(@min,lon)) : min(cellfun(@min,dx)) : max(cellfun(@max,lon));
y_regrid = min(cellfun(@min,lat)) : min(cellfun(@min,dy)) : max(cellfun(@max,lat));
dx = min(cellfun(@min,dx)); dy = min(cellfun(@min,dy));

[xx_regrid,yy_regrid] = meshgrid(x_regrid,y_regrid); %new grid to project data onto

% pre-allocate
vel_regrid = zeros([size(xx_regrid) nframes]);
vstd_regrid = zeros(size(vel_regrid)); mask_regrid = zeros(size(vel_regrid));
compE_regrid = zeros(size(vel_regrid)); compN_regrid = zeros(size(vel_regrid));
compU_regrid = zeros(size(vel_regrid));

for ii = 1:nframes
    
    [xx,yy] = meshgrid(lon{ii},lat{ii});
    vel_regrid(:,:,ii) =  interp2(xx,yy,vel{ii},xx_regrid,yy_regrid);
    vstd_regrid(:,:,ii) =  interp2(xx,yy,vstd{ii},xx_regrid,yy_regrid);
    if par.usemask == 1
        mask_regrid(:,:,ii) =  interp2(xx,yy,mask{ii},xx_regrid,yy_regrid);
    end
    
    % ENU may not be downsampled by licsbas, hence different xx yy grids
    [xx,yy] = meshgrid(lon_comp{ii},lat_comp{ii});
    compE_regrid(:,:,ii) =  interp2(xx,yy,compE{ii},xx_regrid,yy_regrid);
    compN_regrid(:,:,ii) =  interp2(xx,yy,compN{ii},xx_regrid,yy_regrid);
    compU_regrid(:,:,ii) =  interp2(xx,yy,compU{ii},xx_regrid,yy_regrid);
    
    if (mod(ii,round(nframes./10))) == 0
        disp([num2str(round((ii./nframes)*100)) '% completed']);
    end

end

if par.tie2gnss == 1
    [xx_gnss,yy_gnss] = meshgrid(gnss_field.x,gnss_field.y);
    gnss_E = interp2(xx_gnss,yy_gnss,gnss_field.E,xx_regrid,yy_regrid);
    gnss_N = interp2(xx_gnss,yy_gnss,gnss_field.N,xx_regrid,yy_regrid);
end

%% downsample

if par.ds_factor > 0
    
    disp(['Downsampling by a factor of ' num2str(par.ds_factor)])
    
    % get downsampled grid coords
    [~,x_regrid,y_regrid] = downsample_array(vel_regrid(:,:,1),...
        par.ds_factor,par.ds_factor,par.ds_method,x_regrid,y_regrid);
    [xx_regrid,yy_regrid] = meshgrid(x_regrid,y_regrid);
    dx = mean(diff(x_regrid)); dy = mean(diff(y_regrid));
    
    % pre-allocate
    vel_regrid_ds = zeros([length(y_regrid) length(x_regrid) nframes]);
    vstd_regrid_ds = zeros(size(vel_regrid_ds)); mask_regrid_ds = zeros(size(vel_regrid_ds));
    compE_regrid_ds = zeros(size(vel_regrid_ds)); compN_regrid_ds = zeros(size(vel_regrid_ds));
    compU_regrid_ds = zeros(size(vel_regrid_ds));
    
    for ii = 1:nframes
    
        [vel_regrid_ds(:,:,ii),~,~] ...
            = downsample_array(vel_regrid(:,:,ii),par.ds_factor,par.ds_factor,par.ds_method);
        [vstd_regrid_ds(:,:,ii),~,~] ...
            =  downsample_array(vstd_regrid(:,:,ii),par.ds_factor,par.ds_factor,par.ds_method);
        [compE_regrid_ds(:,:,ii),~,~] ...
            =  downsample_array(compE_regrid(:,:,ii),par.ds_factor,par.ds_factor,par.ds_method);
        [compN_regrid_ds(:,:,ii),~,~] ...
            =  downsample_array(compN_regrid(:,:,ii),par.ds_factor,par.ds_factor,par.ds_method);
        [compU_regrid_ds(:,:,ii),~,~] ...
            =  downsample_array(compU_regrid(:,:,ii),par.ds_factor,par.ds_factor,par.ds_method);
        
        if par.usemask == 1
            [mask_regrid_ds(:,:,ii),~,~] ...
                =  downsample_array(mask_regrid(:,:,ii),par.ds_factor,par.ds_factor,par.ds_method);
        end
        
        if (mod(ii,round(nframes./10))) == 0
            disp([num2str(round((ii./nframes)*100)) '% completed']);
        end
        
    end
    
    vel_regrid = vel_regrid_ds; vstd_regrid = vstd_regrid_ds;
    compE_regrid = compE_regrid_ds; compN_regrid = compN_regrid_ds;
    compU_regrid = compU_regrid_ds; mask_regrid = mask_regrid_ds;
    clear vel_regrid_ds vstd_regrid_ds inc_regrid_ds phi_regrid_ds mask_regrid_ds
    
    if par.tie2gnss == 1
        [gnss_E,~,~] = downsample_array(gnss_E,par.ds_factor,par.ds_factor,par.ds_method);
        [gnss_N,~,~] = downsample_array(gnss_N,par.ds_factor,par.ds_factor,par.ds_method);
    end
    
end

%% apply mask

if par.usemask == 1
    
    disp('Applying mask')
    
    mask_regrid(isnan(mask_regrid)) = 0;
    
    % account for previous downsampling
    mask_regrid(mask_regrid>=0.5) = 1; mask_regrid(mask_regrid<0.5) = 0;
    
    vel_regrid(mask_regrid==0) = NaN;
    
    for ii=1:nframes
        vel{ii}(mask{ii}==0) = nan;
    end

end

%% tie to gnss

if par.tie2gnss == 1
    
    % pre-allocate
    gnss_resid_plane = zeros([size(xx_regrid) nframes]);
    
    for ii = 1:nframes
        
        % skip loop if vel is empty (likely because of masking)
        if all(isnan(vel_regrid(:,:,ii)),'all')
            disp([frames{ii} ' vel is empty after masking, skipping referencing'])
            continue
        end
        
        % convert gnss fields to los
        gnss_los = (gnss_E.*compE_regrid(:,:,ii)) + (gnss_N.*compN_regrid(:,:,ii));
        
        % calculate residual
        vel_tmp = vel_regrid(:,:,ii); 
        vel_tmp(vel_tmp>mean(vel_tmp(:),'omitnan')+std(vel_tmp(:),'omitnan')) = nan;
        vel_tmp(vel_tmp<mean(vel_tmp(:),'omitnan')-std(vel_tmp(:),'omitnan')) = nan;
        
        gnss_resid = vel_tmp - gnss_los;
        
        % remove nans
        gnss_xx = xx_regrid(~isnan(gnss_resid));
        gnss_yy = yy_regrid(~isnan(gnss_resid));
        gnss_resid = gnss_resid(~isnan(gnss_resid));
        
        % centre coords
        midx = (max(gnss_xx) + min(gnss_xx))/2;
        midy = (max(gnss_yy) + min(gnss_yy))/2;
        gnss_xx = gnss_xx - midx ;gnss_yy = gnss_yy - midy;
        all_xx = xx_regrid - midx; all_yy = yy_regrid - midy;
          
        % fit plane
        G_resid = [ones(length(gnss_xx),1) gnss_xx gnss_yy gnss_xx.*gnss_yy ...
            gnss_xx.^2 gnss_yy.^2];
        m_resid = (G_resid'*G_resid)^-1*G_resid'*gnss_resid;
        gnss_resid_plane(:,:,ii) = m_resid(1) + m_resid(2).*all_xx + m_resid(3).*all_yy ...
            + m_resid(4).*all_xx.*all_yy + m_resid(5).*all_xx.^2 + m_resid(6).*all_yy.^2;
            
        % remove from insar
        vel_regrid(:,:,ii) = vel_regrid(:,:,ii) - gnss_resid_plane(:,:,ii);
        
    end
    
end

% calculate frame overlaps if requested
if par.frame_overlaps == 1
    disp('Calculating frame overlap statistics')
    frame_overlap_stats(vel_regrid,frames,compU_regrid);
end

%% merge tracks

if par.merge_tracks == 1
    
    % merge frames along tracks by taking mean of overlaps
    [vel_regrid,unique_tracks] = merge_tracks(frames,vel_regrid);
    [vstd_regrid,~] = merge_tracks(frames,vstd_regrid);
    [compE_regrid,~] = merge_tracks(frames,compE_regrid);
    [compN_regrid,~] = merge_tracks(frames,compN_regrid);
    [compU_regrid,~] = merge_tracks(frames,compU_regrid);
    
    nframes = length(unique_tracks);
    
end

%% velocity decomposition

disp('Inverting for E and U')

% pre-al
m_up = nan(size(xx_regrid));
m_east = nan(size(xx_regrid));
m_north = gnss_N;
var_up = nan(size(xx_regrid));
var_east = nan(size(xx_regrid));

npixels = length(xx_regrid(:));

% project gnss into los and remove
for ii = 1:nframes
    gnss_Nlos = gnss_N .* compN_regrid(:,:,ii);
    vel_regrid(:,:,ii) = vel_regrid(:,:,ii) - gnss_Nlos;
end   

% loop through pixels
for jj = 1:size(xx_regrid,1)
    for kk = 1:size(xx_regrid,2)
        
        % make components
        Qd = diag(squeeze(vstd_regrid(jj,kk,:)));
        G = [squeeze(compU_regrid(jj,kk,:)) squeeze(compE_regrid(jj,kk,:))];
        d = squeeze(vel_regrid(jj,kk,:));
        
        % remove invalid pixels
        invalid_pixels = find(isnan(d));
        d(invalid_pixels) = [];
        G(invalid_pixels,:) = [];
        Qd(invalid_pixels,:) = []; Qd(:,invalid_pixels) = [];
        
        % apply cond(G) threshold
        if par.condG_threshold > 0 && cond(G) > par.condG_threshold
            disp('cond(G) threshold exceeded, skipping')
%             m = nan(1,2); Qm = nan(2,2);
            continue
        end
        
        if length(d) < 2 % skip if less than two frames
            continue
        end
        
        % solve
        W = inv(Qd);
        m = (G'*W*G)^-1 * G'*W*d;
        Qm = inv(G'*W*G);
                
        if par.var_threshold > 0 && any(diag(Qm) > par.var_threshold)
            disp('var threshold exceeded, skipping')
%             m = nan(1,2); Qm = nan(2,2);
            continue
        end
        
        
        % save
        m_up(jj,kk) = m(1);
        m_east(jj,kk) = m(2);    
        var_up(jj,kk) = Qm(1,1);
        var_east(jj,kk) = Qm(2,2);
        
    end
    disp([num2str(jj) '/' num2str(size(xx_regrid,1)) ' rows completed'])
end

%% plot output velocities

load('/nfs/a285/homes/eearw/gmt/colourmaps/vik/vik.mat')

lonlim = [min(x_regrid) max(x_regrid)];
latlim = [min(y_regrid) max(y_regrid)];
clim = [-10 10];

f = figure();
f.Position([1 3 4]) = [600 1600 1200];
tiledlayout(2,2,'TileSpacing','compact')

t(1) = nexttile; hold on
plt_data(x_regrid,y_regrid,m_up,lonlim,latlim,clim,'Vertical (mm/yr)',fault_trace,borders)
colormap(t(1),vik)

t(2) = nexttile; hold on
plt_data(x_regrid,y_regrid,m_east,lonlim,latlim,clim,'East (mm/yr)',fault_trace,borders)
colormap(t(2),vik)

t(3) = nexttile; hold on
plt_data(x_regrid,y_regrid,m_north,lonlim,latlim,[],'North (mm/yr)',fault_trace,borders)

t(4) = nexttile; hold on
coverage = sum(~isnan(vel_regrid),3);
plt_data(x_regrid,y_regrid,coverage,lonlim,latlim,[],'Data coverage',fault_trace,borders)

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

%% plotting function ------------------------------------------------------
function plt_data(lon,lat,data,lonlim,latlim,clim,titlestr,fault_trace,borders)

% plot input data
imagesc(lon,lat,data,'AlphaData',~isnan(data));

% plot borders
if ~isempty(borders)
    for ii = 1:length(borders.places)
        plot(borders.lon{ii},borders.lat{ii},'k')
    end
end

% plot fault traces
if ~isempty(fault_trace); plot(fault_trace(:,1),fault_trace(:,2),'r'); end

xlim(lonlim)
ylim(latlim)

colorbar

if ~isempty(clim); caxis(clim); end

title(titlestr)

end