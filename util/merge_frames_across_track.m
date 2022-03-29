function merge_frames_across_track(par,x,y,vel,tracks,compE,compU,vstd)
%=================================================================
% function merge_frames_across_track()
%-----------------------------------------------------------------
% Merge frames across track (assuming they have already been merged
% along-track) by projecting LOS velocities into the average incidence and
% azimuth angle, then removing an offset between adjacent tracks.
%
% These merged velocities are not passed on to the velocity decomposition.
% This is designed as a check of the consistency of the InSAR velocities,
% and so that the structure of the velocity field can be inspected without
% influence from GNSS.
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   tracks: track names (cell array of strings)
%   compE, compU: regridded component vectors (3D arrays)
%   vstd: regridded velocity uncertainties
%   
% Andrew Watson     28-03-2022
%                                                                  
%=================================================================

%% setup

% desired inc and az
av_inc = 40;
av_az_asc = -10;
av_az_desc = -170;

%% get pass dir

% pre-allocate
tracks_asc = cell(0); tracks_desc = cell(0);
tracks_asc_ind = 1:length(tracks); 
tracks_desc_ind = 1:length(tracks); 

% sort tracks by pass direction 
for ii = 1:length(tracks)
    if tracks{ii}(4) == 'A'
        tracks_asc{end+1} = tracks{ii};
        tracks_desc_ind(ii) = 0;
    elseif tracks{ii}(4) == 'D'
        tracks_desc{end+1} = tracks{ii};
        tracks_asc_ind(ii) = 0;
    end
end

tracks_asc_ind(tracks_asc_ind==0) = [];
tracks_desc_ind(tracks_desc_ind==0) = [];

%% calcualte inc and az from components

% get inc from U
inc = acosd(compU);

% then get az from E
az = acosd(compE./sind(inc))-180;

%% project into fixed az and inc
    
vel(:,:,tracks_asc_ind) = vel(:,:,tracks_asc_ind) ...
    .* ( cosd(av_az_asc-az(:,:,tracks_asc_ind)) .* cosd(av_inc-inc(:,:,tracks_asc_ind)) );

vel(:,:,tracks_desc_ind) = vel(:,:,tracks_desc_ind) ...
    .* ( cosd(av_az_desc-az(:,:,tracks_desc_ind)) .* cosd(av_inc-inc(:,:,tracks_desc_ind)) );

%% minimise difference between adjacent tracks
% Assuming tracks are of similar length, order them by min(x).
% Note, this is going to fail if one track is much longer than the other.

% find min x value
min_x = zeros(1,size(vel,3));
for ii = 1:size(vel,3)
    min_x(ii) = x(find(any(vel(:,:,ii)),1));
end

% split into direction
min_x_asc = min_x(tracks_asc_ind);
min_x_desc = min_x(tracks_desc_ind);

% sort
[~,track_order_asc] = sort(min_x_asc);
[~,track_order_desc] = sort(min_x_desc);

% loop through ascending
for ii = 1:length(track_order_asc)-1
    % calculate residual between overlap and remove nans           
    overlap_resid = vel(:,:,tracks_asc_ind(track_order_asc(ii+1))) ...
        - vel(:,:,tracks_asc_ind(track_order_asc(ii)));
    overlap_resid(isnan(overlap_resid)) = [];

    if isempty(overlap_resid)
        disp('No overlap, skipping. Something has gone wrong with track ordering.')
        continue
    end

    % solve linear inverse for offset
    G = ones(length(overlap_resid),1);
    m = (G'*G)^-1*G'*overlap_resid(:);

    % apply offset
    vel(:,:,tracks_asc_ind(track_order_asc(ii+1))) ...
        = vel(:,:,tracks_asc_ind(track_order_asc(ii+1))) - m;    
end

% loop through descending
for ii = 1:length(track_order_desc)-1
    % calculate residual between overlap and remove nans           
    overlap_resid = vel(:,:,tracks_desc_ind(track_order_desc(ii+1))) ...
        - vel(:,:,tracks_desc_ind(track_order_desc(ii)));
    overlap_resid(isnan(overlap_resid)) = [];

    if isempty(overlap_resid)
        disp('No overlap, skipping. Something has gone wrong with track ordering.')
        continue
    end

    % solve linear inverse for offset
    G = ones(length(overlap_resid),1);
    m = (G'*G)^-1*G'*overlap_resid(:);

    % apply offset
    vel(:,:,tracks_desc_ind(track_order_desc(ii+1))) = vel(:,:,tracks_desc_ind(track_order_desc(ii+1))) - m;    
end

%% loop through points

% % coord arrays
% [xx,yy] = meshgrid(x,y);
% 
% % pre-al
% los_av_asc = nan(size(xx)); los_av_desc = nan(size(xx));
% los_vstd_asc = nan(size(xx)); los_vstd_desc = nan(size(xx));
% 
% % number of points in grid
% npixels = length(xx(:));
% 
% % progress report interval
% report_it = round(size(xx,1)/10);
% 
% % loop through pixels
% for jj = 1:size(xx,1)
%     for kk = 1:size(xx,2)
%         
%         % FOR ASCENDING
%         d = squeeze(vel(jj,kk,tracks_asc_ind));
%         G = squeeze(cosd(av_az_asc-az(jj,kk,tracks_asc_ind)) .* cosd(av_inc-inc(jj,kk,tracks_asc_ind)));
%         Qd = diag(squeeze(vstd(jj,kk,tracks_asc_ind)));
%         
%         % remove invalid pixels
%         invalid_pixels = find(isnan(d));
%         d(invalid_pixels) = [];
%         G(invalid_pixels,:) = [];
%         Qd(invalid_pixels,:) = []; Qd(:,invalid_pixels) = [];
%         
%         % solve and save
%         if ~isempty(d)
%             W = inv(Qd);
%             m = (G'*W*G)^-1 * G'*W*d;
%             Qm = inv(G'*W*G);
%             
%             los_av_asc(jj,kk) = m; 
%             los_vstd_asc(jj,kk) = Qm(1,1);
%         end
%                 
%         
%         % FOR DESCENDING
%         d = squeeze(vel(jj,kk,tracks_desc_ind));
%         G = squeeze(cosd(av_az_desc-az(jj,kk,tracks_desc_ind)) .* cosd(av_inc-inc(jj,kk,tracks_desc_ind)));
%         Qd = diag(squeeze(vstd(jj,kk,tracks_desc_ind)));
%         
%         % remove invalid pixels
%         invalid_pixels = find(isnan(d));
%         d(invalid_pixels) = [];
%         G(invalid_pixels,:) = [];
%         Qd(invalid_pixels,:) = []; Qd(:,invalid_pixels) = [];
%         
%         % solve and save
%         if ~isempty(d)
%             W = inv(Qd);
%             m = (G'*W*G)^-1 * G'*W*d;
%             Qm = inv(G'*W*G);
%             
%             los_av_desc(jj,kk) = m; 
%             los_vstd_desc(jj,kk) = Qm(1,1);
%         end       
%         
%     end
%     
%     % report progress
%     if mod(jj,report_it) == 0
%         disp([num2str(jj) '/' num2str(size(xx,1)) ' rows completed'])
%     end
% end

%% take average in overlaps

los_av_asc = mean(vel(:,:,tracks_asc_ind),3,'omitnan');
los_av_desc = mean(vel(:,:,tracks_desc_ind),3,'omitnan');

%% plot

% set plotting parameters
lonlim = [min(x) max(x)];
latlim = [min(y) max(y)];
clim = [-40 40];
load('/nfs/a285/homes/eearw/gmt/colourmaps/vik/vik.mat')

% reload borders for ease
if par.plt_borders == 1
    borders = load(par.borders_file);
else
    borders = [];
end

f = figure();
f.Position([1 3 4]) = [600 1600 600];
tiledlayout(1,2,'TileSpacing','compact')

% plot ascending tracks
t(1) = nexttile; hold on
plt_data(x,y,los_av_asc,lonlim,latlim,clim,'Ascending (mm/yr)',[],borders)
colormap(t(1),vik)

% plot descending tracks
t(2) = nexttile; hold on
plt_data(x,y,los_av_desc,lonlim,latlim,clim,'Descending (mm/yr)',[],borders)
colormap(t(2),vik)

end

