function [vel_regrid,vstd_regrid,mask_regrid] ...
    = apply_masks(par,x_regrid,y_regrid,vel_regrid,vstd_regrid,...
    mask_regrid,asc_frames_ind,desc_frames_ind,fault_trace,borders)
%=================================================================
% function apply_masks()
%-----------------------------------------------------------------
% Apply grid and polygon masks to inputs
%                                                                  
% INPUT:                                                           
%   temp:
% OUTPUT:
%   temp:
%
% Andrew Watson     01-05-2023
%                                                                  
%=================================================================

%% apply masks

[xx_regrid,yy_regrid] = meshgrid(x_regrid,y_regrid);

% tif masks
if par.use_mask == 1 || par.use_mask == 3
    
    disp('Applying tifs masks')
    
    % account for previous downsampling
    mask_regrid(isnan(mask_regrid)) = 0;
    mask_regrid(mask_regrid>=0.5) = 1; mask_regrid(mask_regrid<0.5) = 0;
    vel_regrid(mask_regrid==0) = NaN;
    vstd_regrid(mask_regrid==0) = NaN;
    
    % convert to logical
    mask_regrid = logical(mask_regrid);

end

% poly mask
if par.use_mask == 2 || par.use_mask == 3
    
    disp('Applying poly mask')
    
    % for each polygon
    for ii = 1:length(poly_mask)
        
        % check if within polygon
        [in_poly,~] = inpolygon(xx_regrid(:),yy_regrid(:),poly_mask(ii).X,poly_mask(ii).Y);
        
        % reshape to 3D array
        in_poly = reshape(in_poly,size(xx_regrid));
        in_poly = repmat(in_poly,1,1,nframes);
        
        vel_regrid(in_poly) = nan;
        vstd_regrid(in_poly) = nan;
        
        % report progress
        if (mod(ii,round(length(poly_mask)./10))) == 0
            disp([num2str(round((ii./length(poly_mask))*100)) '% completed']);
        end
        
    end

end

%% plot ascending and descending masks
% Produces two masks that are useful for plotting.
% For overlaps, a given point is considered unmasked if it is unmasked in
% at least one of the masks.

if par.plt_mask_asc_desc == 1
    
    disp('Plotting ascending and descending masks')
    
    % look direction indices
%     asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), frames));
%     desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), frames));
    
    % sum masks
    mask_asc = sum(mask_regrid(:,:,asc_frames_ind),3);
    mask_desc = sum(mask_regrid(:,:,desc_frames_ind),3);
    
    % return to ones and zeros
    mask_asc(mask_asc>=1) = 1;
    mask_desc(mask_desc>=1) = 1;
    mask_asc(mask_asc<1) = 0;
    mask_desc(mask_desc<1) = 0;
    mask_asc(isnan(mask_asc)) = 0;
    mask_desc(isnan(mask_desc)) = 0;
    
    % plot
    lonlim = [min(x_regrid) max(x_regrid)];
    latlim = [min(y_regrid) max(y_regrid)];
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title('Combined masks for ascending and descending')
    
    t(1) = nexttile; hold on
    plt_data(x_regrid,y_regrid,mask_asc,lonlim,latlim,[],'Ascending mask',fault_trace,borders)
    
    t(2) = nexttile; hold on
    plt_data(x_regrid,y_regrid,mask_desc,lonlim,latlim,[],'Descending mask',fault_trace,borders)
    
end

end