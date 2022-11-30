function [lon,lat,dx,dy,lon_comp,lat_comp,vel,vstd,compE,compN,compU,mask,frames,...
    asc_frames_ind,desc_frames_ind,fault_trace,gnss,borders] = load_inputs(par,insarpar)
%=================================================================
% function [] = load_inputs()
%-----------------------------------------------------------------
% Read the input parameter file and load all requested inputs (e.g. frame
% velocities, GNSS velocities).
%                                                                  
% INPUT:                                                           
%   par:  structure containing general parameters
% OUTPUT:                                                          
%   
%
% Andrew Watson     24-11-2022
%                                                                  
%=================================================================

%% setup

% number of velocity maps inputted
nframes = length(insarpar.dir);

% pre-allocate
frames = cell(1,nframes);
lon = cell(1,nframes); lat = cell(size(lon));
lon_comp = cell(size(lon)); lat_comp = cell(size(lon));
dx = cell(size(lon)); dy = cell(size(lon));
vel = cell(size(lon)); vstd = cell(size(lon)); mask = cell(size(lon));
compE = cell(size(lon)); compN = cell(size(lon)); compU = cell(size(lon));

%% load insar

% for each velocity map
for ii = 1:nframes
    
    disp(['Loading ' insarpar.dir{ii}])

    % extract the frame name
    frames(ii) = regexp(insarpar.dir{ii},'\d*[AD]_\d*_\d*','match');

    % load velocities
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vel '*']);
    [lon{ii},lat{ii},vel{ii},dx{ii},dy{ii}] ...
        = read_geotiff([insarpar.dir{ii} namestruct.name],'single');

    % test that vel contains valid pixels, and remove if not
    if sum(~isnan(vel{ii}),'all') == 0
        disp([insarpar.dir{ii} ' is all nans - removing'])
        lon(ii) = []; lat(ii) = []; vel(ii) = []; dx(ii) = []; dy(ii) = [];
        continue
    end

    % load velocity errors
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vstd '*']);
    [~,~,vstd{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name],'single');

    % load East, North, and Up components
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_e '*']);
    [lon_comp{ii},lat_comp{ii},compE{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name],'single');

    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_n '*']);
    [~,~,compN{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name],'single');

    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_u '*']);
    [~,~,compU{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name],'single');
    
    % load mask
    if par.usemask == 1
        namestruct = dir([insarpar.dir{ii} '*' insarpar.id_mask '*']);
        [~,~,mask{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name],'single');
    end
    
end

% get indices of ascending and descending frames
asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), frames));
desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), frames));

%% load gnss

% load gnss vels and check if uncertainties are present (if requested)

switch par.ref2gnss
    case 0
        disp('Not loading GNSS vels - this may break other functions')
        
    case 1
        % GNSS file should be a (at least) six column text file of:
        % lon lat vE vN sE sN (cor)
        
        % check extension
        [~,~,ext] = fileparts(par.gnss_file);
        if strcmp(ext,'.mat')
            error('GNSS stations requested, gnss_file should not end .mat')
        end
        
        gnss = readmatrix(par.gnss_file);
        
    case 2
        
        % check extension
        [~,~,ext] = fileparts(par.gnss_file);
        if ~strcmp(ext,'.mat')
            error('GNSS fields requested, gnss_file should end .mat')
        end
        
        % GNSS file is interpolated velocities in a .mat
        gnss = importdata(par.gnss_file);

        % check for uncertainties
        if par.gnss_uncer == 1 && (~isfield(gnss_field,'sE') || ~isfield(gnss_field,'sN'))
            error('Propagation of GNSS uncertainties requested, but GNSS mat file does not contain sE and sN')
        end
end

%% load plotting files

% fault traces
if par.plt_faults == 1
    fault_trace = single(readmatrix(par.faults_file,'FileType','text'));
else
    fault_trace = [];
end

% borders
if par.plt_borders == 1
    borders = load(par.borders_file);
else
    borders = [];
end

end