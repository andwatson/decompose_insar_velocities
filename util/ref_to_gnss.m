function [vel] = ref_to_gnss(method,xx,yy,vel,compE,compN,gnss_E,gnss_N,frames)
%=================================================================
% function ref_to_gnss()
%-----------------------------------------------------------------
% Tie InSAR velocities into a GNSS refernce frame.
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN, compU: regridded component vectors (3D arrays)
%   vstd: regridded velocity uncertainties
% OUTPUT:    
%   vel: velocities in GNSS reference system
%   
% Andrew Watson     06-06-2022
%                                                                  
%=================================================================

% pre-allocate
nframes = size(vel,3);
gnss_resid_plane = zeros([size(xx) nframes]);

for ii = 1:nframes

    % skip loop if vel is empty (likely because of masking)
    if all(isnan(vel(:,:,ii)),'all')
        disp([frames{ii} ' vel is empty after masking, skipping referencing'])
        continue
    end

    % convert gnss fields to los
    gnss_los = (gnss_E.*compE(:,:,ii)) + (gnss_N.*compN(:,:,ii));

    % calculate residual
    vel_tmp = vel(:,:,ii); 
    vel_tmp(vel_tmp>mean(vel_tmp(:),'omitnan')+std(vel_tmp(:),'omitnan')) = nan;
    vel_tmp(vel_tmp<mean(vel_tmp(:),'omitnan')-std(vel_tmp(:),'omitnan')) = nan;

    gnss_resid = vel_tmp - gnss_los;
    
    % method switch
    switch method
        case 1 % second order polynomial surface
            
            % remove nans
            gnss_xx = xx(~isnan(gnss_resid));
            gnss_yy = yy(~isnan(gnss_resid));
            gnss_resid = gnss_resid(~isnan(gnss_resid));
    
            % centre coords
            midx = (max(gnss_xx) + min(gnss_xx))/2;
            midy = (max(gnss_yy) + min(gnss_yy))/2;
            gnss_xx = gnss_xx - midx ;gnss_yy = gnss_yy - midy;
            all_xx = xx - midx; all_yy = yy - midy;

            % fit plane
            G_resid = [ones(length(gnss_xx),1) gnss_xx gnss_yy gnss_xx.*gnss_yy ...
                gnss_xx.^2 gnss_yy.^2];
            m_resid = (G_resid'*G_resid)^-1*G_resid'*gnss_resid;
            gnss_resid_plane(:,:,ii) = m_resid(1) + m_resid(2).*all_xx + m_resid(3).*all_yy ...
                + m_resid(4).*all_xx.*all_yy + m_resid(5).*all_xx.^2 + m_resid(6).*all_yy.^2;
            
        case 2 % filtering
            
            % filter gnss residual
            gnss_resid_filtered = ndnanfilter(gnss_resid,'rectwin',[121 121]);
            
            % reapply nans
%             gnss_resid_filtered(isnan(gnss_resid)) = nan;
            
            gnss_resid_plane(:,:,ii) = gnss_resid_filtered;
            
    end

    % remove from insar
    vel(:,:,ii) = vel(:,:,ii) - gnss_resid_plane(:,:,ii);

end

end