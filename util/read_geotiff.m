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

if strcmp(georef.RasterInterpretation,'cells')
    disp('Tif is in cell format, converting to postings.')
    dx = georef.CellExtentInLongitude;
    dy = georef.CellExtentInLatitude;
    latlim = [georef.LatitudeLimits(1)+dy/2 georef.LatitudeLimits(2)-dy/2];
    lonlim = [georef.LongitudeLimits(1)+dx/2 georef.LongitudeLimits(2)-dx/2];
    rasterSize = georef.RasterSize;
    georef = georefpostings(latlim,lonlim,rasterSize);
end

% grid spacing
dx = georef.SampleSpacingInLongitude;
dy = georef.SampleSpacingInLatitude;
% try
%     dx = georef.SampleSpacingInLongitude;
%     dy = georef.SampleSpacingInLatitude;
% catch
%     dx = georef.CellExtentInLongitude;
%     dy = georef.CellExtentInLatitude;
% end

% get coords and grid spacing
lon = georef.LongitudeLimits(1) : dx : georef.LongitudeLimits(2);
lat = georef.LatitudeLimits(2) : -dy : georef.LatitudeLimits(1);

end

