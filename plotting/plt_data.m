function plt_data(lon,lat,data,lonlim,latlim,clim,titlestr,fault_trace,borders)
%=================================================================
% function plt_data(lon,lat,data,lonlim,latlim,clim,titlestr,fault_trace,borders)
%-----------------------------------------------------------------
% Plot given array, looping if input is a cell array or 3D.
% Figure must be created outside of function.
%                                                                  
% INPUT:                                                           
%   lon, lat: coordinate vectors
%   data: array or arrays in cell array
%   lonlim, latlim: 1x2 vectors of axes limits
%   clim: 1x2 vector of colour bar limits
%   titlestr: string for title of plot
%   fault_trace: coordinates of fault trace
%   borders: structure of border polygons to add countries to map
%   
% Andrew Watson     23-03-2022
%                                                                  
%=================================================================

hold on

% sparse conversion if required
if issparse(data)
    data = full(data);
    data(data==0) = nan;
end

% plot input data
if iscell(data) % if cell array
    for ii = 1:length(data)
        imagesc(lon{ii},lat{ii},data{ii},'AlphaData',~isnan(data{ii}));
    end
    
elseif ndims(data) == 3 % if 3D array
    for ii = 1:size(data,3)
        imagesc(lon,lat,data(:,:,ii),'AlphaData',~isnan(data(:,:,ii)));
    end
    
else % if 2D array
    imagesc(lon,lat,data,'AlphaData',~isnan(data));
    
end
    
% plot borders
if ~isempty(borders)
    for ii = 1:length(borders.places)
        plot(borders.lon{ii},borders.lat{ii},'k')
    end
end

% plot fault traces
if ~isempty(fault_trace); plot(fault_trace(:,1),fault_trace(:,2),'r'); end

if ~isempty(lonlim); xlim(lonlim); end
if ~isempty(latlim); ylim(latlim); end

colorbar

if ~isempty(clim); caxis(clim); end

title(titlestr)

end