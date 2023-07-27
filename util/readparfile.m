function [par,insarpar] = readparfile(cfgfile)
%=================================================================
% function [parmat] = readparfile(cfgfile)
%-----------------------------------------------------------------
% Function to get parameter values from a config file.
% A description of the available confif options is given in the README.
%                                                                  
% INPUT:                                                           
%   cfgfile: path to parameter text file (e.g. velmap.conf)
% OUTPUT:                                                          
%   par:  structure containing general parameters (paths, processing
%           toggles, plotting and output toggles)
%   insarpar: structure containing InSAR specific options (extensions,
%               number of frames, directories)
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

% gnss fields path
par.gnss_fields_file = getparval(cfgcell,'gnss_fields_file',[]);

% gnss stations path
par.gnss_stations_file = getparval(cfgcell,'gnss_stations_file',[]);

% fault path
par.faults_file = getparval(cfgcell,'faults_file',[]);

% borders path
par.borders_file = getparval(cfgcell,'borders_file',[]);

% external shapefile mask
par.poly_mask_file = getparval(cfgcell,'poly_mask_file',[]);

% output path
par.out_path = getparval(cfgcell,'out_path',[]);

% output prefix
par.out_prefix = getparval(cfgcell,'out_prefix',[]);

%% processing toggles

% regridding coordinates
par.auto_regrid = getparval(cfgcell,'auto_regrid',1);
if par.auto_regrid == 0
    par.crop_post_regrid = getparval(cfgcell,'crop_post_regrid',0);
    
    regrid_xmin = getparval(cfgcell,'regrid_xmin',[]);
    regrid_xmax = getparval(cfgcell,'regrid_xmax',[]);
    regrid_dx = getparval(cfgcell,'regrid_dx',[]);
    regrid_ymin = getparval(cfgcell,'regrid_ymin',[]);
    regrid_ymax = getparval(cfgcell,'regrid_ymax',[]);
    regrid_dy = getparval(cfgcell,'regrid_dy',[]);
    
    par.regrid_x = regrid_xmin:regrid_dx:regrid_xmax;
    par.regrid_y = regrid_ymin:regrid_dy:regrid_ymax;
end

% scale input velocity uncertainties
par.scale_vstd = getparval(cfgcell,'scale_vstd',0);
par.scale_vstd_model = getparval(cfgcell,'scale_vstd_model','sph');

% tie to gnss
par.ref2gnss = getparval(cfgcell,'ref2gnss',0);
par.ref_type = getparval(cfgcell,'ref_type',[]);
par.ref_poly_order = getparval(cfgcell,'ref_poly_order',[]);
par.ref_filter_window_size = getparval(cfgcell,'ref_filter_window_size',[]);
par.ref_station_radius = getparval(cfgcell,'ref_station_radius',0);
par.refmask_min = getparval(cfgcell,'refmask_min',-10);
par.refmask_max = getparval(cfgcell,'refmask_max',10);
par.store_ref_planes = getparval(cfgcell,'store_ref_planes',[]);
par.use_stored_ref_planes = getparval(cfgcell,'use_stored_ref_planes',[]);

% use mask
par.use_mask = getparval(cfgcell,'use_mask',0);

% downsampling
par.ds_factor = getparval(cfgcell,'ds_factor',0);
par.ds_method = getparval(cfgcell,'ds_method','mean',[]);

% marge along-track
par.merge_tracks_along = getparval(cfgcell,'merge_tracks_along',0);
par.merge_tracks_along_func = getparval(cfgcell,'merge_tracks_along_func',0);

% merge across-track
par.merge_tracks_across = getparval(cfgcell,'merge_tracks_across',0);
par.ref_xmin = getparval(cfgcell,'ref_xmin',0)
par.ref_xmax = getparval(cfgcell,'ref_xmax',0)
par.ref_ymin = getparval(cfgcell,'ref_ymin',0)
par.ref_ymax = getparval(cfgcell,'ref_ymax',0)

% reference frame bias
par.plate_motion = getparval(cfgcell,'plate_motion',0);
par.plate_motion_file = getparval(cfgcell,'plate_motion_file',[]);

% gnss uncertainty
par.gnss_uncer = getparval(cfgcell,'gnss_uncer',0);

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

% save grd files
par.save_grd = getparval(cfgcell,'save_grd',0);

% save frames or tracks
par.save_frames = getparval(cfgcell,'save_frames',0);

% save overlaps as text files for plotting histograms
par.save_overlaps = getparval(cfgcell,'save_overlaps',0);

% climits
par.plt_cmin = getparval(cfgcell,'plt_cmin',-10)
par.plt_cmax = getparval(cfgcell,'plt_cmax',10)

% plot faults
par.plt_faults = getparval(cfgcell,'plt_faults',0);

% plot borders
par.plt_borders = getparval(cfgcell,'plt_borders',0);

% plot input vels
par.plt_input_vels = getparval(cfgcell,'plt_input_vels',0);

% plot scaled uncertainties
par.plt_scale_vstd_indv = getparval(cfgcell,'plt_scale_vstd_indv',0);
par.plt_scale_vstd_all = getparval(cfgcell,'plt_scale_vstd_all',0);

% plot merged tracks
par.plt_merge_tracks = getparval(cfgcell,'plt_merge_tracks',0);

% plot plate motion bias corrections
par.plt_plate_motion = getparval(cfgcell,'plt_plate_motion',0);
par.plt_plate_motion_indv = getparval(cfgcell,'plt_plate_motion_indv',0);

% along track merge plotting
par.plt_merge_along_corr = getparval(cfgcell,'plt_merge_along_corr',0);
par.plt_merge_along_resid = getparval(cfgcell,'plt_merge_along_resid',0);

% plot reference to gnss
par.plt_ref_gnss_indv = getparval(cfgcell,'plt_ref_gnss_indv',0);

% plot all referencing surfaces
par.plt_ref_gnss_surfaces = getparval(cfgcell,'plt_ref_gnss_surfaces',0);

% plot gnss in los
par.plt_ref_gnss_los = getparval(cfgcell,'plt_ref_gnss_los',0);

% output gnss in los as grd
par.grd_ref_gnss_los = getparval(cfgcell,'grd_ref_gnss_los',0);

% plot ascending and descending masks
par.plt_mask_asc_desc = getparval(cfgcell,'plt_mask_asc_desc',0);

% plot decomposed velocity uncertainties
par.plt_decomp_uncer = getparval(cfgcell,'plt_decomp_uncer',0);

% plot var and cond(G) threshold masks
par.plt_threshold_masks = getparval(cfgcell,'plt_threshold_masks',0);

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