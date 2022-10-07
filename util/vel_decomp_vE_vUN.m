function [m_east,m_up,var_east,var_up,condG_threshold_mask,var_threshold_mask] ...
    = vel_decomp_vE_vUN(par,vel,vstd,compE,compN,compU,gnss_N,gnss_sN,both_coverage)
%=================================================================
% function vel_decomp()
%-----------------------------------------------------------------
% Decompose InSAR LOS velocities into East and a joint Up-North component.
% The Up-North velocities is the seperately decomposed into Up and North
% using the North GNSS velocitiy.
%                                                                  
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   vel: regridded velocities (3D array)
%   vstd: regridded velocity uncertainties
%   compE, compN, compU: regridded component vectors (3D arrays)
%   gnss_N, gnss_sN: north GNSS velocities and associated uncertainties (2D
%       arrays)
%   both_coverage: logical array true where point has at least one
%       ascending and descending velocity
% OUTPUT:    
%   m_east: decomposed east velocities (2D array)
%   m_up: descomposed vertical velocities
%   var_east: east uncertainty
%   var_up: vertical uncertainty
%   condG_threshold_mask: cond(G) mask (2D array)
%   var_threshold_mask: variance mask
%   
% Andrew Watson     12-09-2022
%                                                                  
%=================================================================

%% setup

% size consts
rowcol = size(vel,[1 2]);
nframes = size(vel,3);

% pre-al
m_UN = nan(rowcol);
m_east = nan(rowcol);
var_UN = nan(rowcol);
var_east = nan(rowcol);
var_threshold_mask = zeros(rowcol);
condG_threshold_mask = zeros(rowcol);

% calculate UN component vector, first by calculating the incidence angle
% and heading. Incidence angle is measured from the vertical, and azimuth
% is measured negatively counterclockwise from north. This is to match the
% definitions in Qi's work.
inc = 90 - asind(compU);
az = acosd(compE./sind(inc))-180;
compUN = sqrt(1 - sind(inc).^2 .* cosd(az).^2);

% number of points in grid
npixels = numel(vel(:,:,1));

%% decompose pixel by pixel

% create loop indexes
[jj,kk] = ndgrid(1:size(vel,1),1:size(vel,2));
jj = jj(:); kk = kk(:);

% reshape array for optimal looping (each pixel becomes a row of a 2D array
% to avoid squeeze within loop).
vel = reshape(vel,[],nframes);
vstd = reshape(vstd,[],nframes);
compUN = reshape(compUN,[],nframes);
compE = reshape(compE,[],nframes);

% progress report interval
report_it = round(size(vel,1)/10);

% loop through pixels
for ii = 1:size(vel,1)
    
    % report progress
    if mod(ii,report_it) == 0
        disp([num2str(ii) '/' num2str(size(vel,1)) ...
            ' (' num2str(round(ii/size(vel,1).*100)) '%) rows completed'])
    end
    
    % skip points without coverage in both look directions
    if both_coverage(jj(ii),kk(ii)) == 0
        continue
    end

    % make components
    Qd = diag(vstd(ii,:));
    G = [compUN(ii,:)' compE(ii,:)'];
    d = vel(ii,:)';

    % remove invalid pixels
    invalid_pixels = find(isnan(d));
    d(invalid_pixels) = [];
    G(invalid_pixels,:) = [];
    Qd(invalid_pixels,:) = []; Qd(:,invalid_pixels) = [];

    % apply cond(G) threshold
    if par.condG_threshold > 0 && cond(G) > par.condG_threshold
        condG_threshold_mask(jj(ii),kk(ii)) = 1;
        m = nan(1,2); Qm = nan(2,2);
        continue
    end

    % solve
    W = inv(Qd);
    m = (G'*W*G)^-1 * G'*W*d;
    Qm = inv(G'*W*G);

    % apply model variance threshold
    if par.var_threshold > 0 && any(diag(Qm) > par.var_threshold)
        var_threshold_mask(jj(ii),kk(ii)) = 1;
        m = nan(1,2); Qm = nan(2,2);
        continue
    end

    % save
    m_UN(jj(ii),kk(ii)) = m(1);
    m_east(jj(ii),kk(ii)) = m(2);    
    var_UN(jj(ii),kk(ii)) = Qm(1,1);
    var_east(jj(ii),kk(ii)) = Qm(2,2);
    
end

% report number of points removed.
disp([num2str(sum(condG_threshold_mask,'all')) '/' num2str(npixels) ...
    ' (' num2str(round(sum(condG_threshold_mask,'all')/npixels*100),2) ...
    '%) points were masked by the cond(G) threshold.'])
disp([num2str(sum(var_threshold_mask,'all')) '/' num2str(npixels) ...
    ' (' num2str(round(sum(var_threshold_mask,'all')/npixels*100),2) ...
    '%) points were masked by the model variance threshold.'])

%% decompose vUN into vU and vN

% estimate vU for all frames/tracks
UN2U = (sqrt(1 - sind(inc).^2 .* cosd(az).^2)) ./ cosd(inc);
N2U = (sind(az).*sind(inc)) ./ cosd(inc);
m_up = m_UN.*UN2U - gnss_N.*N2U;

% same for the uncertainty
UN2U_var = (1 - sind(inc).^2 .* cosd(az).^2) ./ cosd(inc);
N2U_var = ((sind(az).^2).*(sind(inc).^2)) ./ (cosd(inc).^2);
var_up = sqrt((var_UN.^2 .* UN2U_var) + (gnss_sN .* N2U_var));

% take the weighted mean
m_up = sum(m_up.*(1./var_up),3,'omitnan') ./ sum((1./var_up),3,'omitnan');
var_up = sum(var_up.^2,3,'omitnan') ./ sum(var_up,3,'omitnan');


