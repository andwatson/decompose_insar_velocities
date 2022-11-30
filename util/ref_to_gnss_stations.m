function [vel] ...
    = ref_to_gnss_stations(par,cpt,xx,yy,vel,compE,compN,gnss,asc_frames_ind,desc_frames_ind)
%=================================================================
% function ref_to_gnss_stations()
%-----------------------------------------------------------------
% Tie InSAR velocities into a GNSS refernce frame.
% Use individual GNSS station velocities, as opposed to the interpolated
% fields used in ref_to_gnss_fields (original method).
%
% Main steps are:
%   - project GNSS into LOS for each frame/track
%   - calculate residual between InSAR and GNSS
%   - fit a polynomial to the residual
%   - subtract this polynomial from the InSAR
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   cpt: structure containing colour palettes
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN, compU: regridded component vectors (3D arrays)
%   gnss: matrix of gnss station locations, velocities, and uncertainties
%   asc_frames_ind, desc_frames_ind: indices for ascending and descending
%       frames/tracks
% OUTPUT:    
%   vel: velocities in GNSS reference system
%   
% Andrew Watson     30-11-2022
%                                                                  
%=================================================================

nframes = size(vel,3);
nstations = size(gnss,1);

%% check for sufficient gnss coverage
% There must be at least as many gnss stations within the velocity fields
% for each frame/track as there are parameters in the polynomial, otherwise
% the inverse problem cannot be solved.

% number of stations required
n_stations_required = par.ref_poly_order .* 3;

for ii = 1:nframes
    
    % get coords for valid pixels
    xx_check = xx(~isnan(vel(:,:,ii)));
    yy_check = yy(~isnan(vel(:,:,ii)));
    
    % distance between all stations and every valid point
    xx_check = repmat(xx_check,1,nstations);
    yy_check = repmat(yy_check,1,nstations);
    
    dists = sqrt( (xx_check - gnss(:,1)').^2 + (yy_check - gnss(:,2)').^2 );
    
    % number of stations with at least one vel within radius
    n_stations_within = sum(sum(dists <= par.ref_station_radius) >= 1);
    
    if n_stations_within < n_stations_required
        error(['Vel ' num2str(ii) ' contains insufficient gnss stations to perform referening'])
    end
    
end

%% perform referencing

% first, get indices of coords within radius of each station
within_radius = cell(1,nstations);
for jj = 1:nstations
    dists = sqrt( (xx - gnss(jj,1)').^2 + (yy - gnss(jj,2)').^2 );
    within_radius{jj} = find(dists <= par.ref_station_radius);
end

% remove any stations with no vels within radius
gnss(cellfun(@isempty, within_radius),:) = [];
within_radius(cellfun(@isempty, within_radius)) = [];
nstations = length(within_radius);

% pre-al
gnss_resid = zeros(nstations,nframes);
gnss_resid_poly = nan(size(vel));

% reference
for ii = 1:nframes
    
%     gnss_los = zeros(nstations,1);
    
    % project to LOS and calculate residual
    for jj = 1:nstations
        
        [within_row,within_col] = ind2sub(size(xx),within_radius{jj});
        
        gnss_los = (gnss(jj,3) .* mean(compE(within_row,within_col,ii),'all','omitnan')) ...
            + (gnss(jj,4) .* mean(compN(within_row,within_col,ii),'all','omitnan'));
        
        gnss_resid(jj,ii) = mean(vel(within_row,within_col,ii),'all','omitnan') - gnss_los;
                
    end
    
    % remove nans
    resid = gnss_resid(~isnan(gnss_resid(:,ii)),ii);
    x = gnss(~isnan(gnss_resid(:,ii)),1);
    y = gnss(~isnan(gnss_resid(:,ii)),2);
    
    % fit polynomial to residuals
    switch par.ref_poly_order
        case 1
            G_resid = [ones(length(x),1) x y];
            m_resid = (G_resid'*G_resid)^-1*G_resid'*resid; % ADD UNCER
            gnss_resid_poly(:,:,ii) = m_resid(1) + m_resid(2).*xx + m_resid(3).*yy;
            
    end
    
    
end


% mask resid with vel (just for plotting)
vel_mask = single(~isnan(vel(:,:,ii)));
vel_mask(vel_mask==0) = nan;
gnss_resid_plane(:,:,ii) = gnss_resid_poly(:,:,ii) .* vel_mask;

% remove from insar
vel(:,:,ii) = vel(:,:,ii) - gnss_resid_poly(:,:,ii);

end