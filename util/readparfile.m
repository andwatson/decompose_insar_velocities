function [par,insarpar] = readparfile(cfgfile)
%=================================================================
% function [parmat] = readparfile(cfgfile)
%-----------------------------------------------------------------
% Function to get parameter values from a configure file
%                                                                  
% INPUT:                                                           
%   cfgfile: path to parameter text file (e.g. velmap.conf)
% OUTPUT:                                                          
%   par:  structure containing general parameters
%   
% Adapted from the velmap function of the same name.
%
% Andrew Watson     20-08-2021
%                                                                  
% NOTE: use '%' for comments in config file, and ': ' to seperate names and
% values (e.g. inv_e:   1)
%=================================================================

%% open config file

% test that config file exists
if ~isfile(cfgfile)
    disp('Config file not found, continuing with default values.')
end

% load the config file as a cell array
cfgcell = readcell(cfgfile,'FileType','text','CommentStyle','%','Delimiter',':');

%% paths

% gnss path
par.gnss_file = getparval(cfgcell,'gnss_file',[]);

% fault path
par.faults_file = getparval(cfgcell,'faults_file',[]);

% borders path
par.borders_file = getparval(cfgcell,'borders_file',[]);

% output path
par.out_path = getparval(cfgcell,'out_path',[]);

% output prefix
par.out_prefix = getparval(cfgcell,'out_prefix',[]);

%% processing toggles

% tie to gnss
par.tie2gnss = getparval(cfgcell,'tie2gnss',0);

% use mask
par.usemask = getparval(cfgcell,'usemask',0);

% downsampling
par.ds_factor = getparval(cfgcell,'ds_factor',0);
par.ds_method = getparval(cfgcell,'ds_method','mean',[]);

% marge along-track
par.merge_tracks_along = getparval(cfgcell,'merge_tracks_along',0);
par.merge_tracks_along_func = getparval(cfgcell,'merge_tracks_along_func',0);

% merge across-track
par.merge_tracks_across = getparval(cfgcell,'merge_tracks_across',0);

% reference frame bias
par.plate_motion = getparval(cfgcell,'plate_motion',0);
par.plate_motion_file = getparval(cfgcell,'plate_motion_file',[]);

% decomposition method
par.decomp_method = getparval(cfgcell,'decomp_method',0);

% threshold on cond(G) in inversion
par.condG_threshold = getparval(cfgcell,'condG_threshold',0);

% threshold on model variance in inversion
par.var_threshold = getparval(cfgcell,'var_threshold',0);

% calculate frame overlap statistics
par.frame_overlaps = getparval(cfgcell,'frame_overlaps',0);

%% plotting and output toggles

% save geotiffs
par.save_geotif = getparval(cfgcell,'save_geotif',0);

% plot faults
par.plt_faults = getparval(cfgcell,'plt_faults',0);

% plot borders
par.plt_borders = getparval(cfgcell,'plt_borders',0);

% plot input vels
par.plt_input_vels = getparval(cfgcell,'plt_input_vels',0);

% plot merged tracks
par.plt_merge_tracks = getparval(cfgcell,'plt_merge_tracks',0);

%% insar

% file identifiers
insarpar.id_vel = getparval(cfgcell,'id_vel','vel');
insarpar.id_vstd = getparval(cfgcell,'id_vstd','vstd');
insarpar.id_e = getparval(cfgcell,'id_e','E');
insarpar.id_n = getparval(cfgcell,'id_n','N');
insarpar.id_u = getparval(cfgcell,'id_u','U');
insarpar.id_mask = getparval(cfgcell,'id_mask','mask');
insarpar.id_par = getparval(cfgcell,'id_par','par');

% insar data
insarpar.ninsarfile = sum(strcmp(cfgcell(:,1),'framedir'));
for ii = 1:insarpar.ninsarfile
    insarpar.dir{ii} = getparval(cfgcell,'framedir',[],ii);
end


%% getparval ==============================================================
% cfgcell: n-by-2 cell array containing parameters
% parstring: name of desired parameter
% defval: default val if parameter is not found
% indn: index to use when multiple matches are found
function [val] = getparval(cfgcell,parstring,defval,indn)
    
    if nargin ~= 4
        indn = 1;
    end
    
    try
%         val = cfgcell{strcmp(cfgcell(:,1), parstring),2};
        indx = find(strcmp(cfgcell(:,1), parstring));
        val = cfgcell{indx(indn),2};
    catch
        val = defval;
    end
end

end