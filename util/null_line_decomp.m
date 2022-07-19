function [m_perp1,m_perp2,var_east,var_up,condG_threshold_mask,var_threshold_mask] ...
    = null_line_decomp(par,vel,vstd,compE,compN,compU,both_coverage,asc_frames_ind,desc_frames_ind)
%=================================================================
% function null_line_decomp()
%-----------------------------------------------------------------
% Decompose InSAR LOS velocities into the two vectors orthogonal to the
% null line.
%                                                                  
% INPUT:                                                           
%   par: input parameter structure
%   vel: 3D array of velocities
%   vstd: uncertainties on velocities
%   compE, compN, compU: LOS component vectors
%   both_coverage: mask for where at least two look directions are present
%   asc_frames_ind: indices of ascending frames
%   desc_frames_ind: indicies of descending frames
% OUTPUT:    
%   m_perp1: decomp vel in first perp direction
%   m_perp2: decomp vel in second perp direction
%   var_east: variance for perp one 
%   var_up: variance for perp two
%   condG_threshold_mask: mask based on cond(G) threshold
%   var_threshold_mask: mask based on variance threshold
%   
% Andrew Watson     13-06-2022
%                                                                  
%=================================================================

%% setup

% size consts
rowcol = size(vel,[1 2]);
nframes = size(vel,3);

% pre-al
m_perp1 = nan(rowcol);
m_perp2 = nan(rowcol);
var_up = nan(rowcol);
var_east = nan(rowcol);
var_threshold_mask = zeros(rowcol);
condG_threshold_mask = zeros(rowcol);

% number of points in grid
npixels = numel(vel(:,:,1));

%% find null line and perp vectors
% assume a single null line for the entire area

% average of components for asc and desc
comp_asc = [mean(compE(:,:,asc_frames_ind),'all','omitnan')
    mean(compN(:,:,asc_frames_ind),'all','omitnan')
    mean(compU(:,:,asc_frames_ind),'all','omitnan')];

comp_desc = [mean(compE(:,:,desc_frames_ind),'all','omitnan')
    mean(compN(:,:,desc_frames_ind),'all','omitnan')
    mean(compU(:,:,desc_frames_ind),'all','omitnan')];

% components of perpendicular vectors to null line
[comp_perp1,comp_perp2] = null_line(comp_asc,comp_desc);

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
compN = reshape(compN,[],nframes);
comp_perp1 = comp_perp1(:);
comp_perp2 = comp_perp2(:);

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
%     Qd = diag(vstd(ii,:));
    G = [compE(ii,:)' compN(ii,:)' compU(ii,:)'];
    d = vel(ii,:)';

    % remove invalid pixels
    invalid_pixels = find(isnan(d));
    d(invalid_pixels) = [];
    G(invalid_pixels,:) = [];
%     Qd(invalid_pixels,:) = []; Qd(:,invalid_pixels) = [];

    % convert to null line perps
    G_perp = zeros(size(G,1),2);
    for mm = 1:size(G,1)
        G_perp(mm,1) = dot(G(mm,:)',comp_perp1);
        G_perp(mm,2) = dot(G(mm,:)',comp_perp2);
    end

    % solve
%     W = inv(Qd);
%     m = (G'*W*G)^-1 * G'*W*d;
    m = G_perp\d;
%     Qm = inv(G'*W*G);

%     % apply model variance threshold
%     if par.var_threshold > 0 && any(diag(Qm) > par.var_threshold)
%         var_threshold_mask(jj(ii),kk(ii)) = 1;
%         m = nan(1,2); Qm = nan(2,2);
%         continue
%     end

    % save
    m_perp1(jj(ii),kk(ii)) = m(1);
    m_perp2(jj(ii),kk(ii)) = m(2);    
%     var_up(jj(ii),kk(ii)) = Qm(1,1);
%     var_east(jj(ii),kk(ii)) = Qm(2,2);
    
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