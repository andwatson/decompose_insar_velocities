function plt_asc_desc_cells(par,lon,lat,data,mask,asc_frames_ind,desc_frames_ind,cpt,clim,borders)
%=================================================================
% function preview_input_vels(cfgfile)
%-----------------------------------------------------------------
% Plot input velocities.
%                                                                  
% INPUT:                                                           
%   lon,lat: cell arrays containing coordinate vectors for each frame
%   data: cell array containing 2D data arrays for each frame (vel, vstd)
%   mask: cell array containing the 2D mask arrays for each frame
%   asc_frames_ind,desc_frames_ind: indices for the ascending and
%                                    descending frames
%   nframes: number of input velocity fields
%   cpt: colour palette for plotting
%   borders: structure containing polygons defining country borders 
%
% Andrew Watson     28-04-2023
%                                                                  
%=================================================================

% set plotting parameters
lonlim = [min(cellfun(@min,lon)) max(cellfun(@max,lon))]; 
latlim = [min(cellfun(@min,lat)) max(cellfun(@max,lat))];

% temporarily apply mask for plotting
if par.use_mask == 1
    for ii = 1:length(data)
        data{ii}(mask{ii}==0) = nan;
        % data{ii} = data{ii} - median(data{ii},'all','omitnan');
    end
end

f = figure();
f.Position([1 3 4]) = [600 1600 600];
t = tiledlayout(1,2,'TileSpacing','compact');
title(t,'Input velocities')

% plot ascending tracks
t(1) = nexttile; hold on
plt_data(lon(asc_frames_ind),lat(asc_frames_ind),data(asc_frames_ind),...
    lonlim,latlim,clim,'Ascending (mm/yr)',[],borders)
colormap(t(1),cpt)

% plot descending tracks
t(2) = nexttile; hold on
plt_data(lon(desc_frames_ind),lat(desc_frames_ind),data(desc_frames_ind),...
    lonlim,latlim,clim,'Descending (mm/yr)',[],borders)
colormap(t(2),cpt)

clear vel_tmp


end