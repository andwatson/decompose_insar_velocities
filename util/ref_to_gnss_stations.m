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
%   gnss: structure containing "stations" = matrix of gnss station locations, velocities, and uncertainties
%   asc_frames_ind, desc_frames_ind: indices for ascending and descending
%       frames/tracks
% OUTPUT:    
%   vel: velocities in GNSS reference system
%   
% Andrew Watson     30-11-2022
%                                                                  
%=================================================================

nframes = size(vel,3);
nstations = size(gnss.stations,1);

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
    
    dists = sqrt( (xx_check - gnss.stations(:,1)').^2 + (yy_check - gnss.stations(:,2)').^2 );
    
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
    dists = sqrt( (xx - gnss.stations(jj,1)').^2 + (yy - gnss.stations(jj,2)').^2 );
    within_radius{jj} = find(dists <= par.ref_station_radius);
end

% remove any stations with no vels within radius
gnss.stations(cellfun(@isempty, within_radius),:) = [];
within_radius(cellfun(@isempty, within_radius)) = [];
nstations = length(within_radius);

% calculate los uncertainties
gnss_los_sd = sqrt(gnss.stations(:,5).^2 + gnss.stations(:,6).^2);

% pre-al
gnss_resid = zeros(nstations,nframes);
gnss_resid_poly = nan(size(vel));

% reference
for ii = 1:nframes
    
    % project to LOS and calculate residual
    for jj = 1:nstations
        
        [within_row,within_col] = ind2sub(size(xx),within_radius{jj});
        
        gnss_los = (gnss.stations(jj,3) .* mean(compE(within_row,within_col,ii),'all','omitnan')) ...
            + (gnss.stations(jj,4) .* mean(compN(within_row,within_col,ii),'all','omitnan'));
       
        gnss_resid(jj,ii) = mean(vel(within_row,within_col,ii),'all','omitnan') - gnss_los;
                
    end
    
    % remove nans
    resid = gnss_resid(~isnan(gnss_resid(:,ii)),ii);
    resid_x = gnss.stations(~isnan(gnss_resid(:,ii)),1);
    resid_y = gnss.stations(~isnan(gnss_resid(:,ii)),2);
    resid_W = inv(diag(gnss_los_sd(~isnan(gnss_resid(:,ii)))));
    
    % fit polynomial to residuals
    switch par.ref_poly_order
        case 1
            G_resid = [ones(length(resid_x),1) resid_x resid_y];
            m_resid = (G_resid'*resid_W*G_resid)^-1*G_resid'*resid_W*resid;
            gnss_resid_poly(:,:,ii) = m_resid(1) + m_resid(2).*xx + m_resid(3).*yy;
            
        case 2
            G_resid = [ones(length(resid_x),1) resid_x resid_y resid_x.*resid_y resid_x.^2 resid_y.^2];
            m_resid = (G_resid'*resid_W*G_resid)^-1*G_resid'*resid_W*resid;
            gnss_resid_poly(:,:,ii) = m_resid(1) + m_resid(2).*xx + m_resid(3).*yy ...
                + m_resid(4).*xx.*yy + m_resid(5).*xx.^2 + m_resid(6).*yy.^2;
            
    end
    
    % for plotting
    if par.plt_ref_gnss_indv == 1
        vel_orig = vel(:,:,ii);
    end
    
    % mask resid with vel (just for plotting)
    vel_mask = single(~isnan(vel(:,:,ii)));
    vel_mask(vel_mask==0) = nan;
    gnss_resid_poly(:,:,ii) = gnss_resid_poly(:,:,ii) .* vel_mask;

    % remove from insar
    vel(:,:,ii) = vel(:,:,ii) - gnss_resid_poly(:,:,ii);
    
    % optional plotting
    if par.plt_ref_gnss_indv == 1

        % limits
        clim = [-10 10];
        x = xx(1,:); y = yy(:,1);
        [~,x_ind,y_ind] = crop_nans(vel(:,:,ii),x,y);
        lonlim = x([x_ind(1) x_ind(end)]); latlim = y([y_ind(1) y_ind(end)]);
        
        f = figure();
        f.Position([1 3 4]) = [600 1600 600];
        tiledlayout(1,3,'TileSpacing','compact');

        nexttile; hold on
        imagesc(x,y,vel_orig,'AlphaData',~isnan(vel_orig)); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Original InSAR vel')

        nexttile; hold on
        imagesc(x,y,gnss_resid_poly(:,:,ii),'AlphaData',~isnan(gnss_resid_poly(:,:,ii))); axis xy
        scatter(resid_x,resid_y,40,resid,'Filled','MarkerEdgeColor','k')
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Referencing residual')

        nexttile; hold on
        imagesc(x,y,vel(:,:,ii),'AlphaData',~isnan(vel(:,:,ii))); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Referenced InSAR')
        
    end
    
    % report progress
    disp([num2str(ii) '/' num2str(nframes) ' complete'])
    
end

%% plotting referencing surfaces

if par.plt_ref_gnss_surfaces == 1
    % plot referencing functions
    x = xx(1,:); y = yy(:,1);
    lonlim = [min(x) max(x)];
    latlim = [min(y) max(y)];
    clim = [-10 10];

    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'Referencing surfaces')

    % ascending tracks
    t(1) = nexttile; hold on
    plt_data(x,y,gnss_resid_poly(:,:,asc_frames_ind),lonlim,latlim,clim,'Ascending (mm/yr)',[],[])
    colormap(t(1),cpt.vik)

    % descending tracks
    t(2) = nexttile; hold on
    plt_data(x,y,gnss_resid_poly(:,:,desc_frames_ind),lonlim,latlim,clim,'Descending (mm/yr)',[],[])
    colormap(t(2),cpt.vik)
end

end