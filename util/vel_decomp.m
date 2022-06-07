function [m_east,m_up,var_east,var_up,condG_threshold_mask,var_threshold_mask] ...
    = vel_decomp(par,vel,vstd,compE,compN,compU,gnss_N,gnss_sN,both_coverage)
%=================================================================
% function vel_decomp()
%-----------------------------------------------------------------
% Decompose InSAR LOS velocities into East and Vertical components.
%                                                                  
% INPUT:                                                           
%   par: 
% OUTPUT:    
%   track_vel: 
%   
% Andrew Watson     07-06-2022
%                                                                  
%=================================================================

%% setup

% size consts
rowcol = size(vel,[1 2]);
nframes = size(vel,3);

% pre-al
m_up = nan(rowcol);
m_east = nan(rowcol);
var_up = nan(rowcol);
var_east = nan(rowcol);
var_threshold_mask = zeros(rowcol);
condG_threshold_mask = zeros(rowcol);

% number of points in grid
npixels = numel(vel(:,:,1));

% project gnss into los and remove
for ii = 1:nframes
    gnss_Nlos = gnss_N .* compN(:,:,ii);
    vel(:,:,ii) = vel(:,:,ii) - gnss_Nlos;
    
    % propagate error on N
    vstd(:,:,ii) = sqrt(vstd(:,:,ii).^2 + gnss_sN.^2);
end   

%% decompose pixel by pixel

% create loop indexes
[jj,kk] = ndgrid(1:size(vel,1),1:size(vel,2));
jj = jj(:); kk = kk(:);

% reshape array for optimal looping (each pixel becomes a row of a 2D array
% to avoid squeeze within loop).
vel = reshape(vel,[],nframes);
vstd = reshape(vstd,[],nframes);
compU = reshape(compU,[],nframes);
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
    Qd = diag(vstd(ii,:)');
    G = [compU(ii,:)' compE(ii,:)'];
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
    m_up(jj(ii),kk(ii)) = m(1);
    m_east(jj(ii),kk(ii)) = m(2);    
    var_up(jj(ii),kk(ii)) = Qm(1,1);
    var_east(jj(ii),kk(ii)) = Qm(2,2);
    
end

% report number of points removed.
disp([num2str(sum(condG_threshold_mask,'all')) '/' num2str(npixels) ...
    ' (' num2str(round(sum(condG_threshold_mask,'all')/npixels*100),2) ...
    '%) points were masked by the cond(G) threshold.'])
disp([num2str(sum(var_threshold_mask,'all')) '/' num2str(npixels) ...
    ' (' num2str(round(sum(var_threshold_mask,'all')/npixels*100),2) ...
    '%) points were masked by the model variance threshold.'])

% reshape back to 3D arrays
% vel = reshape(vel,[rowcol nframes]);
% vstd = reshape(vstd,[rowcol nframes]);
% compU = reshape(compU,[rowcol nframes]);
% compE = reshape(compE,[rowcol nframes]);