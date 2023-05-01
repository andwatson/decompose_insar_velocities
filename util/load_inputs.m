function [lon,lat,dx,dy,lon_comp,lat_comp,vel,vstd,compE,compN,compU,mask,poly_mask,frames,...
    asc_frames_ind,desc_frames_ind,fault_trace,gnss,borders] = load_inputs(par,insarpar)
%=================================================================
% function [] = load_inputs()
%-----------------------------------------------------------------
% Read the input parameter file and load all requested inputs (e.g. frame
% velocities, GNSS velocities).
%                                                                  
% INPUT:                                                           
%   par:  structure containing general parameters
%   insarpar: structure containing insar-specific parameters
% OUTPUT:                                                          
%   lon,lat: cell arrays containing coordinate vectors for each frame
%   dx,dy: cell arrays containing input grid spacing for each frame
%   lon_comp,lat_comp: cell arrays containing coordinate vectors for the
%                       look vector component maps for each frame (issue 
%                       with an older version of lics where the vectors
%                       were on the wrong grid, may not be obsolete)
%   vel: cell array containing the 2D velocity arrays for each frame
%   vstd: cell array containing the 2D uncertainty arrays for each frame
%   compE,compN,compU: cell arrays containing the 2D look vector component
%                       arrays for each frame
%   mask: cell array containing the 2D mask arrays for each frame
%   poly_mask: structure array containing the polygons from the poly_mask
%               shape file
%   frames: cell array of frame name strings
%   asc_frames_ind,desc_frames_ind: indices for the ascending and
%                                    descending frames
%   fault_trace: nx2 coordinate array defining fault traces for plotting
%   gnss: structure that optionally combines both the interpolated GNSS
%           fields and the GNSS stations
%   borders: structure containing polygons defining country borders
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
    
    % load mask tifs for each vel
    if par.use_mask == 1 || par.use_mask == 3          
        namestruct = dir([insarpar.dir{ii} '*' insarpar.id_mask '*']);
        [~,~,mask{ii},~,~] = read_geotiff([insarpar.dir{ii} namestruct.name],'single');        
    end
    
    % load shapefile mask
    if par.use_mask == 2 || par.use_mask == 3        
        poly_mask = shaperead(par.poly_mask_file);
    else
        poly_mask = [];      
    end
    
end

% get indices of ascending and descending frames
asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), frames));
desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), frames));

%% load gnss

% select which to load.
% we only need the stations vels if they're used for referencing.
% we need the field velocities if they are used for referencing or in the
% decomp.

% no referencing, and assuming N is zero
if par.ref2gnss == 0 && par.decomp_method == 3
    load_stations = 0; load_fields = 0;
    
% ref to stations, assume N is zero
elseif par.ref2gnss == 1 && par.decomp_method == 3
    load_stations = 1; load_fields = 0;
    
% ref to stations, any other decomp method
elseif par.ref2gnss == 1 && ismember(par.decomp_method,0:2)
    load_stations = 1; load_fields = 1;
    
% ref to fields, any decomp
elseif par.ref2gnss == 2
    load_stations = 0; load_fields = 1;
    
% anything else
else
    load_stations = 0; load_fields = 0;
end

% load neither
if load_stations == 0 && load_fields == 0
    disp('Loading neither GNSS stations or GNSS fields.')
    gnss = [];
end

% load fields
if load_fields == 1
    
    disp('Loading interpolated GNSS fields')
    
    % check extension
    [~,~,ext] = fileparts(par.gnss_fields_file);
    if ~strcmp(ext,'.mat')
        error('GNSS fields requested, gnss_file should end .mat')
    end

    % GNSS file is interpolated velocities in a .mat
    gnss = importdata(par.gnss_fields_file);

    % check for uncertainties
    if par.gnss_uncer == 1 && (~isfield(gnss,'sE') || ~isfield(gnss,'sN'))
        error('Propagation of GNSS uncertainties requested, but GNSS mat file does not contain sE and sN')
    end
end

% load stations
if load_stations == 1
    % GNSS file should be a (at least) six column text file of:
    % lon lat vE vN sE sN (cor)
    
    disp('Loading GNSS station velocities')

    % check extension
    [~,~,ext] = fileparts(par.gnss_stations_file);
    if strcmp(ext,'.mat')
        error('GNSS stations requested, gnss_file should not end .mat')
    end

    gnss.stations = readmatrix(par.gnss_stations_file);
end

% have to load fields first so that the structures combine easily.

%% load misc

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