function [data_crop,x_ind,y_ind,x_crop,y_crop] = crop_nans(data,x,y)
%=================================================================
% function crop_nans()
%-----------------------------------------------------------------
% Crop all-nan rows and columns from input array by indexing based on first
% and last all-nan rows and columns.
%                                                                  
% INPUT:                                                           
%   data:
%   x:
%   y:
% OUTPUT:    
%   data_crop:
%   x_crop:
%   y_crop:
%   
% Andrew Watson     21-06-2022
%                                                                  
%=================================================================

% find rows and columns that contain any non-nan value
val_cols = find(any(~isnan(data),1));
val_rows = find(any(~isnan(data),2));

x_ind = val_cols(1):val_cols(end);
y_ind = val_rows(1):val_rows(end);

% crop data
data_crop = data(y_ind,x_ind);

% optinally crop x and y coords
if nargin == 3
    x_crop = x(x_ind);
    y_crop = y(y_ind);
end