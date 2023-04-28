function [lon,lat,vel,vstd,compE,compN,compU,mask] ...
    = downsamp(par,nframes,lon,lat,vel,vstd,compE,compN,compU,mask)
%=================================================================
% function [] = downsamp()
%-----------------------------------------------------------------
% Downsample vel, vstd, mask, and the unit vector component cell arrays.
%                                                                  
% INPUT:                                                           
%   par:  structure containing general parameters
% OUTPUT:                                                          
%   
%
% Andrew Watson     30-11-2022
%                                                                  
%=================================================================

for ii = 1:nframes

    [vel{ii},lon{ii},lat{ii}] ...
        = downsample_array(vel{ii},par.ds_factor,par.ds_factor,par.ds_method,lon{ii},lat{ii});        
    [vstd{ii},~,~] = downsample_array(vstd{ii},par.ds_factor,par.ds_factor,par.ds_method);

    [compE{ii},lon{ii},lat{ii}] ...
        = downsample_array(compE{ii},par.ds_factor,par.ds_factor,par.ds_method,lon{ii},lat{ii});
    [compN{ii},~,~] = downsample_array(compN{ii},par.ds_factor,par.ds_factor,par.ds_method);
    [compU{ii},~,~] = downsample_array(compU{ii},par.ds_factor,par.ds_factor,par.ds_method);

    if par.use_mask == 1
        [mask{ii},~,~] ...
            =  downsample_array(mask{ii},par.ds_factor,par.ds_factor,par.ds_method);
    end

    if (mod(ii,round(nframes./10))) == 0
        disp([num2str(round((ii./nframes)*100)) '% completed']);
    end

end

end