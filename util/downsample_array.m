function [out_array,out_x,out_y] = downsample_array(in_array,nlook_r,nlook_c,ds_fun,x,y)
%% DOWNSAMPLE_ARRAY
% Provide matrix, vectors of x and y coords, and multilook factor.
% Downsamples matrix by given amount using 'blockproc' from the image
% processing toolbox.
% In the case where the size of the array does not perfectly divide by the
% downsampling factor, blockproc still uses those values to make a new
% cell. That means that one edge will get extrapolated outwards from the
% real data (as opposed to the licsbas multilook method that cuts off those
% extra elements.
%
% downsample function options : 'mean', 'median'
%
% Andrew Watson     02-05-2020

%% DOWNSAMPLE

% downsampling function
switch ds_fun
    case 'mean'
        fun = @(block_struct) mean(block_struct.data(:),'omitnan');
        
    case 'median'
        fun = @(block_struct) median(block_struct.data(:),'omitnan');
end

% downsample main array
out_array = blockproc(in_array,[nlook_r,nlook_c],fun);

if nargin == 6

    %force to be horizontal
    y = reshape(y,1,[]); x = reshape(x,1,[]);

    %pad with extra values if needed
    if mod(length(y),nlook_r) > 0
        diff_tmp = diff(y);
        y = [y y(end)+cumsum(diff_tmp(1:nlook_r-mod(length(y),nlook_r)))];
    end
    if mod(length(x),nlook_c) > 0
        diff_tmp = diff(x);
        x = [x x(end)+cumsum(diff_tmp(1:nlook_c-mod(length(x),nlook_c)))];
    end

    out_y = mean(reshape(y,nlook_c,[]),1);
    out_x = mean(reshape(x,nlook_c,[]),1);

else
    out_y = nan; out_x = nan;
end

end

