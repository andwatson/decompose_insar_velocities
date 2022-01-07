function [lon,lat,data,dx,dy] = read_geotiff(tif_file)
%% read_geotiff.m
% Read geotif file and output data, coord vectors, and grid spacing.
% Note that lat counts downwards
%
% Andrew Watson     24-08-2021

% open geotiff
[data,georef] = readgeoraster(tif_file);

% convert to double
data = double(data);

% get coords and grid spacing
lon = georef.LongitudeLimits(1) : georef.SampleSpacingInLongitude : georef.LongitudeLimits(2);
lat = georef.LatitudeLimits(2) : -georef.SampleSpacingInLatitude : georef.LatitudeLimits(1);
dx = georef.SampleSpacingInLongitude;
dy = georef.SampleSpacingInLatitude;

end

