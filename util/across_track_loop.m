function [across_track_overlaps] = across_track_loop(combins,tracks,unique_tracks,vel)
%=================================================================
% function [] = frame_overlap_stats()
%-----------------------------------------------------------------
% Function to get parameter values from a configure file
%                                                                  
% INPUT:                                                           
%   cfgfile: path to parameter text file (e.g. velmap.conf)
% OUTPUT:                                                          
%   par:  structure containing general parameters
%   
% Adapted from the velmap function of the same name.
%
% Andrew Watson     20-08-2021
%                                                                  
% NOTE: use '%' for comments in config file, and ': ' to seperate names and
% values (e.g. inv_e:   1)
%=================================================================

% pre-allocate
across_track_overlaps = [];

% loop through each combination of overlaps
for ii = 1:size(combins,1)
    
    % get frame inds for the two tracks
    track1_indn = find(cellfun(@(x) strncmp(unique_tracks{combins(ii,1)},x,length(unique_tracks{combins(ii,1)})), tracks));
    track2_indn = find(cellfun(@(x) strncmp(unique_tracks{combins(ii,2)},x,length(unique_tracks{combins(ii,1)})), tracks));
    
    % all combinations of frames in these two tracks
    [c1,c2] = meshgrid(track1_indn,track2_indn);
    frame_combs = [c1(:) c2(:)];
    
    % difference betweeen all frames in tracks
    frame_diff = nan([size(vel,[1 2]) size(frame_combs,1)]);
    for jj = 1:size(frame_combs)
        frame_diff(:,:,jj) = vel(:,:,frame_combs(jj,1)) - vel(:,:,frame_combs(jj,2));
    end
    
    % append to output
    frame_diff = frame_diff(:);
    frame_diff(isnan(frame_diff)) = [];
    across_track_overlaps = [across_track_overlaps; frame_diff];
    
end

end

