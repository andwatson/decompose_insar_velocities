function frame_overlap_stats(vel,frames,compU)
%=================================================================
% function frame_overlap_stats(vel,frames,compU)
%-----------------------------------------------------------------
% Calculate the differences between overlapping frames both along-track and
% across-track, for ascending and descending frames. Plots histograms of
% the overlap differences for all frames, including a gaussian fit to
% estimate the mean and standard deviation.
%
% Based on frame overlap work in Weiss et al. (2020).
%                                                                  
% INPUT:                                                           
%   vel: matrix of regridded velocities
%   frame: cell array of frame names
%   compU: Up component vectors for each vel
%   
% Andrew Watson     04-03-2022
%                                                                  
%=================================================================

% remove nans from vels
vel(isnan(vel)) = 0;

%% get unique tracks, split by pass direction

% get tracks from frame ids, removing duplicates
tracks = cellfun(@(x) x(1:4), frames, 'UniformOutput', false);
unique_tracks = unique(tracks);

% sort tracks by pass direction 
unique_tracks_asc = cell(0); unique_tracks_desc = cell(0);
for ii = 1:length(unique_tracks)
    if unique_tracks{ii}(4) == 'A'
        unique_tracks_asc{end+1} = unique_tracks{ii};
    elseif unique_tracks{ii}(4) == 'D'
        unique_tracks_desc{end+1} = unique_tracks{ii};
    end
end

%% along-track overlaps (both asc and desc)

disp('Calculating along-track overlaps')

% pre-allocate
along_track_asc = nan(size(vel(:,:,1)));
along_track_desc = nan(size(vel(:,:,1)));
track_diffs = cell(1,length(unique_tracks));

% loop through each track
for ii = 1:length(unique_tracks)
    
    % get inds for frames on that track
    track_ind = cellfun(@(x) strncmp(unique_tracks{ii},x,length(unique_tracks{ii})), tracks);
    track_indn = find(track_ind);
    
    % find the overlaps
    overlaps = sum( (vel(:,:,track_ind)~=0 & ~isnan(vel(:,:,track_ind)) ),3);
    
    % subtract each frame
    frame_sub = vel(:,:,track_indn(1));
    for jj = 2:length(track_indn)
        frame_sub = frame_sub - vel(:,:,track_indn(jj));
    end
    
    % combine into asc and desc
    if unique_tracks{ii}(4) == 'A'
        along_track_asc(overlaps==2) = frame_sub(overlaps==2);
    elseif unique_tracks{ii}(4) == 'D'
        along_track_desc(overlaps==2) = frame_sub(overlaps==2);
    end
    
    % store just diffs for each track
    track_diffs{ii} = frame_sub(overlaps==2);
    
end

%% prep across-track overlaps

disp('Calculating across-track overlaps')

% calc inc from u
inc = asind(compU);
inc(inc==0) = nan;

% assume that all motion is horizontal, using the mean incidence angle of
% each frame
vel_horz = vel./cosd(inc);
for ii = 1:length(tracks)
    vel_horz(:,:,ii) = vel_horz(:,:,ii) .* mean(cosd(inc(:,:,ii)),'all','omitnan');
end

% get all combinations of tracks (split into asc and desc)
combins_asc = nchoosek(1:length(unique_tracks_asc),2);
combins_desc = nchoosek(1:length(unique_tracks_desc),2);

% calculate across track overlaps
[across_track_asc] = across_track_loop(combins_asc,tracks,unique_tracks_asc,vel);
[across_track_desc] = across_track_loop(combins_desc,tracks,unique_tracks_desc,vel);

%% plot histograms
% Histogram each set and fit a gaussian function to estimate mean and SD.

disp('Plotting frame overlap histograms')

xmin = -20; xmax = 20; xint = 0.1;
% nhist = 41; hist_edges = linspace(xmin,xmax,nhist);
nhist = 50;

figure()
tiledlayout(2,2,'TileSpacing','compact')

% along-track asc
nexttile; hold on
h1 = histogram(along_track_asc(:),nhist);
h1_mids = h1.BinEdges(1:end-1)+(diff(h1.BinEdges)./2);
f1 = fit(h1_mids',h1.Values','gauss1');
plot(xmin:xint:xmax,f1(xmin:xint:xmax),'r')
legend(['mean = ' num2str(f1.b1)],['std = ' num2str(f1.c1)])
% xlim([xmin xmax])
title('Along-track asc')

% along-track desc
nexttile; hold on
h2 = histogram(along_track_desc(:),nhist);
h2_mids = h2.BinEdges(1:end-1)+(diff(h2.BinEdges)./2);
f2 = fit(h2_mids',h2.Values','gauss1');
plot(xmin:xint:xmax,f2(xmin:xint:xmax),'r')
legend(['mean = ' num2str(f2.b1)],['std = ' num2str(f2.c1)])
% xlim([xmin xmax])
title('Along-track desc')

% across-track asc
nexttile; hold on
h3 = histogram(across_track_asc(:),nhist);
h3_mids = h3.BinEdges(1:end-1)+(diff(h3.BinEdges)./2);
f3 = fit(h3_mids',h3.Values','gauss1');
plot(xmin:xint:xmax,f3(xmin:xint:xmax),'r')
legend(['mean = ' num2str(f3.b1)],['std = ' num2str(f3.c1)])
% xlim([xmin xmax])
title('Across-track asc')

% across-track desc
nexttile; hold on
h4 = histogram(across_track_desc(:),nhist);
h4_mids = h4.BinEdges(1:end-1)+(diff(h4.BinEdges)./2);
f4 = fit(h4_mids',h4.Values','gauss1');
plot(xmin:xint:xmax,f4(xmin:xint:xmax),'r')
legend(['mean = ' num2str(f4.b1)],['std = ' num2str(f4.c1)])
% xlim([xmin xmax])
title('Across-track desc')

end

