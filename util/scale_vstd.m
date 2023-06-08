function [vstd_scaled,rms_misfit] = scale_vstd(par,x,y,vstd,frame)
%=================================================================
% function vstd_scale()
%-----------------------------------------------------------------
% Scale InSAR velocity uncertainties using semivariogram to mitigate the
% impact of the reference pixel.
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   x, y: coord vectors
%   vstd: velocity uncertainties (2D array)
%
% OUTPUT:    
%   vstd_scale: scaled uncertainty (array, same size as vstd)
%   rms_misfit: rms misfit between the model and vstd
%   
% Andrew Watson     19-10-2022
%                                                                  
%=================================================================

%% setup

% convert to UTM
dczone = utmzone(mean(x),mean(y));
utmstruct = defaultm('utm'); 
utmstruct.zone = dczone;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

[xx,yy] = meshgrid(x,y);

[xx_utm,yy_utm] = mfwdtran(utmstruct,yy(:),xx(:));
coords_utm = [xx_utm(:) yy_utm(:)]./1000;

% hardcode number of bins
nbins = 100;

% identfy reference pixel by assuming that min value occurs there
[~,min_ind] = min(vstd(:));
% [ref_row,ref_col] = ind2sub(size(vstd),min_ind);

% calculate distance between reference pixel and all others pixels
% [xx,yy] = meshgrid(x,y);
% ref_dists = sqrt( (xx - xx(ref_row,ref_col)).^2 + (yy - yy(ref_row,ref_col)).^2 );
ref_dists = sqrt( (coords_utm(:,1) - coords_utm(min_ind,1)).^2 ...
    + (coords_utm(:,2) - coords_utm(min_ind,2)).^2 );

% reshape and mask distances with vstd
ref_dists = reshape(ref_dists,size(xx));
ref_dists(isnan(vstd)) = nan;

% calcualte bin medians and SD (hardcoded 100 bins), starting from 0
max_dist = max(ref_dists(:));
bin_edges = linspace(0,max_dist,nbins+1);
bin_mids = bin_edges(1:end-1) + (diff(bin_edges)./2);

bin_medians = zeros(1,nbins);
bin_stds = zeros(1,nbins);

for ii = 1:nbins
    bin_medians(ii) = median(vstd( ref_dists>bin_edges(ii) & ref_dists<=bin_edges(ii+1) ),'all','omitnan');
    bin_stds(ii) = std(vstd( ref_dists>bin_edges(ii) & ref_dists<=bin_edges(ii+1) ),'omitnan');
end

%% fit spherical model

% generate weights
W = 1 ./ (bin_stds + bin_mids./max_dist);

% define objective function to minimise
switch par.scale_vstd_model
    case 'sph'
        obj_fun = @(m) sum(((spherical_model(m,bin_mids)-bin_medians).^2).*W,'omitnan');
    case 'exp'
        obj_fun = @(m) sum(((exponential_model(m,bin_mids)-bin_medians).^2).*W,'omitnan');
end

% starting values
% m0 = [1 0.5 0.1];
m0 = [20 0.5 0.1];

% define limits, main one being that range can't be greater than the max
% profile distance
lower_lim = [0 0 0];
upper_lim = [max_dist./3 inf inf];

% create options for fminsearch
options = optimset('MaxFunEvals',1000000);

% solve
% [m_fit,~,~,~] = fminsearch(obj_fun,m0,options);
[m_fit,~,~,~] = fminsearchbnd(obj_fun,m0,lower_lim,upper_lim,options);

%% scale

switch par.scale_vstd_model
    case 'sph'
        vstd_scaled = vstd .* (m_fit(2)./spherical_model(m_fit,ref_dists));
        bin_fit = spherical_model(m_fit,bin_mids);
        point_fit = spherical_model(m_fit,ref_dists(:));
    case 'exp'
        vstd_scaled = vstd .* (m_fit(2)./exponential_model(m_fit,ref_dists));
        bin_fit = exponential_model(m_fit,bin_mids);
        point_fit = exponential_model(m_fit,ref_dists(:));
end

% calculate rms
rms_misfit = rms(vstd(:)-point_fit,'omitnan');

%% plot original and scaled

if par.plt_scale_vstd_indv == 1

    clim = [min(vstd(:)) max(vstd(:))];

    f = figure();
    f.Position([1 3 4]) = [100 1400 1000];

    tiledlayout(2,2,'TileSpacing','compact');

    % original vario
    nexttile; hold on
    scatter(ref_dists(:),vstd(:),0.1,'MarkerFaceColor',[110 199 38]./255,'MarkerEdgeColor','k');
    plot(bin_mids,bin_medians,'Color',[198 55 188]./255,'LineWidth',2)
    plot(bin_mids,bin_medians+bin_stds,'Color',[198 55 188]./255,'LineWidth',1)
    plot(bin_mids,bin_medians-bin_stds,'Color',[198 55 188]./255,'LineWidth',1)
    plot(bin_mids,bin_fit,'Color','k','LineWidth',2)
    xlabel('Distance')
    ylabel('vstd difference')
    title('Original uncertainty profile')

    % scaled
    nexttile;
    scatter(ref_dists(:),vstd_scaled(:),0.1,'MarkerFaceColor',[110 199 38]./255,'MarkerEdgeColor','k')
    xlabel('Distance')
    ylabel('vstd difference')
    title('Scaled uncertainty profile')

    % original vstd
    nexttile
    imagesc(x,y,vstd,'AlphaData',~isnan(vstd))
    axis xy
    colorbar; caxis(clim);
    title('Original uncertainties')

    % scaled vstd
    nexttile
    imagesc(x,y,vstd_scaled,'AlphaData',~isnan(vstd_scaled))
    axis xy
    colorbar; caxis(clim);
    title('Scaled uncertainties')

end

%% save output

% outdir = '/nfs/a285/homes/eearw/gmt/thesis/chp4/scaled_uncer_all/data/';
% 
% grdwrite2(x,y,vstd,[outdir frame '_pre.grd'])
% grdwrite2(x,y,vstd_scaled,[outdir frame '_post.grd'])
% 
% writematrix([ref_dists(:) vstd(:)],[outdir frame '_points_pre.txt'])
% writematrix([ref_dists(:) vstd_scaled(:)],[outdir frame '_points_post.txt'])
% writematrix([bin_mids(:) bin_fit(:)],[outdir frame '_fit.txt'])
% writematrix(m_fit,[outdir frame '_fit_par.txt'])

%% spherical model

    function out = spherical_model(m,dists)
        out = m(3) + m(2) .* ( ((3.*dists)./(2.*m(1))) - (dists.^3./(2.*m(1).^3)) );
        out(dists>m(1)) = m(3) + m(2);
    end

%% exponential model

    function out = exponential_model(m,dists)
        out = m(3) + m(2).*(1-exp(-dists./m(1)));        
    end

end