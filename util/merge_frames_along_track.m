function [track_vel,track_compE,track_compN,track_compU,track_vstd,unique_tracks] ...
    = merge_frames_along_track(par,x,y,vel,frames,compE,compN,compU,vstd)
%=================================================================
% function merge_frames_along_track()
%-----------------------------------------------------------------
% Merge frames along each track using overlaps, solving for either a static
% offset (merge_tracks_along==1) or a 1st order plane
% (merge_tracks_along==2). 
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN, compU: regridded component vectors (3D arrays)
%   vstd: regridded velocity uncertainties
% OUTPUT:    
%   track_vel: merged velocities for each track (3D array)
%   track_compE, track_compN, track_compU: average component vectors
%   track_vstd: merged velocity uncertainties
%   unique_tracks: track names (cell array of strings)
%   
% Andrew Watson     04-03-2022
%                                                                  
%=================================================================

% pause loop to display frames for each track before and after shift
plt_before_after = 0;

%% get unique tracks, split by pass direction

% get tracks from frame ids, removing duplicates
tracks = cellfun(@(x) x(1:4), frames, 'UniformOutput', false);
unique_tracks = unique(tracks);

% pre-allocate
unique_tracks_asc = cell(0); unique_tracks_desc = cell(0);
unique_tracks_asc_ind = 1:length(unique_tracks); 
unique_tracks_desc_ind = 1:length(unique_tracks); 

% sort tracks by pass direction 
for ii = 1:length(unique_tracks)
    if unique_tracks{ii}(4) == 'A'
        unique_tracks_asc{end+1} = unique_tracks{ii};
        unique_tracks_desc_ind(ii) = 0;
    elseif unique_tracks{ii}(4) == 'D'
        unique_tracks_desc{end+1} = unique_tracks{ii};
        unique_tracks_asc_ind(ii) = 0;
    end
end

unique_tracks_asc_ind(unique_tracks_asc_ind==0) = [];
unique_tracks_desc_ind(unique_tracks_desc_ind==0) = [];

%% merge frames along-track

% pre-allocate
track_vel = nan([size(vel,[1 2]) length(unique_tracks)]);
track_compE = nan(size(track_vel)); track_compN = nan(size(track_vel));
track_compU = nan(size(track_vel)); track_vstd = nan(size(track_vel));

% coords
[xx,yy] = meshgrid(x,y);

% loop through each track
for ii = 1:length(unique_tracks)
    
    % get inds for frames on that track
    track_ind = find(cellfun(@(x) strncmp(unique_tracks{ii},x,length(unique_tracks{ii})), tracks));
        
    switch par.merge_tracks_along
        case 1 % static offset
            
            % loop through adjcent pairings along-track
            for jj = 1:length(track_ind)-1
                
                % calculate residual between overlap and remove nans           
                overlap_resid = vel(:,:,track_ind(jj+1)) - vel(:,:,track_ind(jj));
                overlap_resid(isnan(overlap_resid)) = [];
                
                if isempty(overlap_resid)
                    disp('No overlap, skipping')
                    continue
                end
                
                % solve linear inverse for offset
                G = ones(length(overlap_resid),1);
                m = (G'*G)^-1*G'*overlap_resid(:);
                
                % apply offset
                vel(:,:,track_ind(jj+1)) = vel(:,:,track_ind(jj+1)) - m;
                
            end
            
            
        case 2 % 1st order plane
            
            % loop through adjcent pairings along-track
            for jj = 1:length(track_ind)-1
                
                % calculate residual between overlap and remove nans           
                overlap_resid = vel(:,:,track_ind(jj)) - vel(:,:,track_ind(jj+1));
                x_overlap = xx(~isnan(overlap_resid)); y_overlap = yy(~isnan(overlap_resid));
                overlap_resid(isnan(overlap_resid)) = [];
                
                if isempty(overlap_resid)
                    disp('No overlap, skipping')
                    continue
                end
                
                G = [ones(length(x_overlap),1) x_overlap y_overlap];
                m = (G'*G)^-1*G'*overlap_resid';
                overlap_plane = m(1) + m(2).*xx + m(3).*yy;
                
                % apply offset
                vel(:,:,track_ind(jj)) = vel(:,:,track_ind(jj)) - overlap_plane;
                
            end
            
    end
    
    % take mean of overlap and store new vel
    track_vel(:,:,ii) = mean(vel(:,:,track_ind),3,'omitnan');
    
    % merge component vectors
    track_compE(:,:,ii) = mean(compE(:,:,track_ind),3,'omitnan');
    track_compN(:,:,ii) = mean(compN(:,:,track_ind),3,'omitnan');
    track_compU(:,:,ii) = mean(compU(:,:,track_ind),3,'omitnan');
    
    % merge pixel uncertainties
    track_vstd(:,:,ii) = mean(vstd(:,:,track_ind),3,'omitnan');
    
    
    if plt_before_after == 1
        % plot original frames
        f1 = figure();
        tiledlayout(1,length(track_ind),'TileSpacing','compact')
        for kk = 1:length(track_ind)
            nexttile
            col_ind = find(any(~isnan(vel(:,:,track_ind(kk))),1));
            row_ind = find(any(~isnan(vel(:,:,track_ind(kk))),2));        
            imagesc(x,y,vel(:,:,track_ind(kk)),'AlphaData',~isnan(vel(:,:,track_ind(kk))))
            xlim([x(col_ind(1)) x(col_ind(end))])
            ylim([y(row_ind(1)) y(row_ind(end))])
            caxis([-10 10])
            colorbar
            axis xy
        end

        % plot merged result
        f2 = figure();
        col_ind = find(any(~isnan(track_vel(:,:,ii)),1));
        row_ind = find(any(~isnan(track_vel(:,:,ii)),2));
        imagesc(x,y,track_vel(:,:,ii),'AlphaData',~isnan(track_vel(:,:,ii)));
        xlim([x(col_ind(1)) x(col_ind(end))])
        ylim([y(row_ind(1)) y(row_ind(end))])
        caxis([-10 10])
        colorbar
        axis xy

        waitfor(f1); waitfor(f2)
    end
    
    % report progress
    disp([num2str(ii) '/' num2str(length(unique_tracks)) ' complete'])
    
end


%% plot merged tracks

if par.plt_merge_tracks == 1
    
    % set plotting parameters
    lonlim = [min(x) max(x)];
    latlim = [min(y) max(y)];
    clim = [-10 10];
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
    plt_data(x,y,track_vel(:,:,unique_tracks_asc_ind),lonlim,latlim,clim,'Ascending (mm/yr)',[],borders)
    colormap(t(1),vik)
    
    % plot descending tracks
    t(2) = nexttile; hold on
    plt_data(x,y,track_vel(:,:,unique_tracks_desc_ind),lonlim,latlim,clim,'Descending (mm/yr)',[],borders)
    colormap(t(2),vik)

end

end


% % fit plane
% if ref_surface == 1 % (c x y)
%     G_resid = [ones(length(gnss_xx),1) gnss_xx gnss_yy];
%     m_resid = (G_resid'*G_resid)^-1*G_resid'*gnss_resid;
%     gnss_resid_plane(:,:,ii) = m_resid(1) + m_resid(2).*all_xx + m_resid(3).*all_yy;
% 
% elseif ref_surface == 2 % (c x y xy x^2 y^2)
%     G_resid = [ones(length(gnss_xx),1) gnss_xx gnss_yy gnss_xx.*gnss_yy ...
%         gnss_xx.^2 gnss_yy.^2];
%     m_resid = (G_resid'*G_resid)^-1*G_resid'*gnss_resid;
%     gnss_resid_plane(:,:,ii) = m_resid(1) + m_resid(2).*all_xx + m_resid(3).*all_yy ...
%         + m_resid(4).*all_xx.*all_yy + m_resid(5).*all_xx.^2 + m_resid(6).*all_yy.^2;


%     % get combinations of frames
%     pairings = nchoosek(track_ind,2);
%     overlap_ind = zeros(length(track_ind)-1,2);
%     
%     % find overlaps by testing every combination for non-nan values
%     n = 1;
%     for kk = 1:size(pairings,1)
%         if any(~isnan(vel(:,:,pairings(kk,1)) - vel(:,:,pairings(kk,2))),'all')
%             overlap_ind(n,:) = pairings(kk,:);
%             n = n + 1;
%         end
%     end    