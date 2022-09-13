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

%% get unique tracks, split by pass direction

% get indices of ascending and descending frames (for if not merging)
asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), frames));
desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), frames));

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
if par.plt_merge_along_resid == 1; overlaps = cell(1,size(vel,3)); n_ov = 1; end

% coords
[xx,yy] = meshgrid(x,y);

% loop through each track
for ii = 1:length(unique_tracks)
    
    disp(['Merging ' unique_tracks{ii}])
    
    % get inds for frames on that track
    track_ind = find(cellfun(@(x) strncmp(unique_tracks{ii},x,length(unique_tracks{ii})), tracks));
        
    switch par.merge_tracks_along_func
        case 0 % static offset
            
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
                
                if par.plt_merge_along_resid == 1
                    overlaps{n_ov} = vel(:,:,track_ind(jj+1)) - vel(:,:,track_ind(jj));
                    n_ov = n_ov + 1;
                end
                
            end
            
            
        case 1 % 1st order plane
            
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
                
                if par.plt_merge_along_resid == 1
                    overlaps{n_ov} = vel(:,:,track_ind(jj+1)) - vel(:,:,track_ind(jj));
                    n_ov = n_ov + 1;
                end
                
            end
            
    end
    
    % plot residuals between overlaps
    if par.plt_merge_along_resid == 1
        % pre-al
        overlaps = overlaps(~cellfun('isempty',overlaps));
        overlap_stats = zeros(length(overlaps),2);
        
        % loop through overlaps
        for nn = 1:length(overlaps)
            figure(); tiledlayout(2,1,'TileSpacing','compact')
            [cropped_overlap,~,~,x_crop,y_crop] = crop_nans(overlaps{nn},x,y);
            nexttile(); imagesc(x_crop,y_crop,cropped_overlap,'AlphaData',~isnan(cropped_overlap)); 
            colorbar; axis xy
            nexttile(); histogram(overlaps{nn});
            title(['Mean = ' num2str(round(mean(cropped_overlap(:),'omitnan'),2)) ...
                ', SD = ' num2str(round(std(cropped_overlap(:),'omitnan'),2)) ...
                'Median = ' num2str(round(median(cropped_overlap(:),'omitnan'),2))])
            overlap_stats(nn,:) = [mean(cropped_overlap(:),'omitnan') std(cropped_overlap(:),'omitnan')];
        end
    end
    
    % merge frames in track into single velocity field
    if par.merge_tracks_along == 2
        % take mean of overlap and store new vel
        track_vel(:,:,ii) = mean(vel(:,:,track_ind),3,'omitnan');

        % merge component vectors
        track_compE(:,:,ii) = mean(compE(:,:,track_ind),3,'omitnan');
        track_compN(:,:,ii) = mean(compN(:,:,track_ind),3,'omitnan');
        track_compU(:,:,ii) = mean(compU(:,:,track_ind),3,'omitnan');

        % merge pixel uncertainties
        track_vstd(:,:,ii) = mean(vstd(:,:,track_ind),3,'omitnan');
    end
    
    % plot original velocities and merged result
    if par.plt_merge_along_corr == 1
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

% if not merging overlaps, then duplicate variables for output
if par.merge_tracks_along == 1
    track_vel = vel; track_vstd = vstd;
    track_compE = compE; track_compN = compN; track_compU = compU;
    unique_tracks_asc_ind = asc_frames_ind; unique_tracks_desc_ind = desc_frames_ind;
end

%% plot merged tracks

if par.plt_merge_tracks == 1
    
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
 