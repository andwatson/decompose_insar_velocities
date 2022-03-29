function tif2grd(in_file,out_file)
%=================================================================
% function tif2grd(in_file,out_file)
%-----------------------------------------------------------------
% Convert a geotiff file to a gmt grd file.                                                                
%   
% Andrew Watson     07-03-2022
%                                                                  
%=================================================================

% automate out_file name if not given
if nargin == 1
    [pathstr, name, ~] = fileparts(in_file);
    out_file = [pathstr '/' name '.grd'];
end

% load input tif
[lon,lat,data,~,~] = read_geotiff(in_file);

% write output grd
grdwrite2(lon,lat,data,out_file);

end

