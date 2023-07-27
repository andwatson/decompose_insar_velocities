function [data_deramp] = deramp_fixed_heading(x,y,data,theta)
%=================================================================
% function deramp(x,y,data,method)
%-----------------------------------------------------------------
% Deramp using a 1st order polynomial plane that tilts about a given
% heading angle. This was originally written to remove azimuth ramps
% observed in the frame velocities. We now understand these to be plate
% motion effects, and suggest using "plate_motion_bias.m" to remove these
% instead.
%                                                                  
% INPUT:                                                           
%   x, y: coordinates
%   data: array or column vector of data
%   theta: heading angle
% OUTPUT:    
%   data_deramp: deramped data
%   
% Andrew Watson     02-06-2022
%                                                                  
%================================================================

% plot toggle
plt = 0;

% duplicate (partially for plotting)
data_orig = data;
data_deramp = data;

% grid coords if needed
if numel(x) ~= numel(data)
    [xx,yy] = meshgrid(x,y);
else
    xx = x; yy = y;
end

% strip nans
xx(isnan(data)) = [];
yy(isnan(data)) = [];
data(isnan(data)) = [];

% solve as linear inverse
G = [ones(numel(xx),1) cosd(theta).*xx(:)+sind(theta).*yy(:)];
m = G\data';

% evaluate and remove
[xx,yy] = meshgrid(x,y);
ramp = m(1)+m(2).*(cosd(theta).*xx+sind(theta).*yy);
data_deramp = data_deramp - ramp;

%% plot

if plt == 1

    % set plotting parameters
    lonlim = [min(x) max(x)];
    latlim = [min(y) max(y)];
    clim = [par.plt_cmin par.plt_cmax];
    load('vik.mat')

    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    tiledlayout(1,3,'TileSpacing','compact')

    % plot original
    t(1) = nexttile; hold on
    plt_data(x,y,data_orig,lonlim,latlim,clim,'Original',[],[])
    colormap(t(1),vik)

    % plot deramp
    t(2) = nexttile; hold on
    plt_data(x,y,ramp,lonlim,latlim,clim,'Ramp',[],[])
    contour(xx,yy,ramp,'k')
    colormap(t(2),vik)

    % plot deramped
    t(3) = nexttile; hold on
    plt_data(x,y,data_deramp,lonlim,latlim,clim,'Deramp',[],[])
    colormap(t(3),vik)

end

