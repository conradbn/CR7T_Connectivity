%% CALL SURFACE PLOTTING FUNCTION
% Benjamin Neal Conrad, Vanderbilt University, November 2018

clear all; close all
% Add subfolders to path  
addpath(genpath(pwd))

%% PLOT VALUED DATA TO SURFACE EXAMPLE
% Read in data (vector of values corresponding to each ROI in the label/surface files)
data = dlmread('test_data_valued.txt');
data = data(:,1); % Just plot the first row (volume) of data

% Base colormap
c = jet(256);
% RdBu colormap (requires cbrewer package)
% RdBu = cbrewer('div', 'RdBu', 256);
cmap_base = c;

% Color range (min val, max val)
colorrange = [-1.5,1.5];

% Surface gifti file
fs_dir = 'fsaverage_BrainnetomeAtlas/164k/' ;
surface = [fs_dir 'fsaverage.R.very_inflated.164k.surf.gii'];

% Label gifti file
label = [fs_dir 'fsaverage.R.BN_Atlas.164k.label.gii'];

% Output directory
out_dir = pwd;

% Output file name
out_names = {'Surf_Valued_lat_r_90_0';...
             'Surf_Valued_med_r_270_0';...
             'Surf_Valued_sup_r_0_90'};

% Viewing angle
view_angs = [90,0;...
            270,0;...
            0,90];
        
% ----------- RUN FUNCTION --------------
for ii = 1:numel(out_names)
    out_name = out_names{ii};
    view_ang = view_angs(ii,:);
    plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);
end


%-------------------------------------------------------------------------------------------------
%% PLOT INTEGER DATA TO SURFACE EXAMPLE
cmap_base = ones(256,3);
% Set module/network integer value, followed by the desired RGB
% Note the module labels aren't used in the present script
module_id_rgb = {'DMN',2,1,0,0;... % DMN, red
                 'FPN',6,1,1,0;... % FPN, yellow
                 'SMN',10,0,1,1;... % Sensorimotor, cyan
                 'Sal',12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
                 'STS',16,1,0.7,1;... % STS, dark gray
                 'CaudV',20,1,0.5,0;... % ventral Caudate, orange
                 'Vis',21,0,0,1;... % Visual, blue
                 'Hipp',22,0,0.5,0.5;... % Hippocampus, teal
                 'DAN',24,0,1,0;... % Dorsal Attention (ITG/MTG),green
                 'Prec',27,0.5,0.5,1;... % Precuneus, pale blue
                 'BG',34,0.3,0.1,0;... % Basal ganglia, dark brown
                 'CaudD',35,.8,.33,0;... % Caudate, burnt orange
                 'Thal',36,0.5,0.25,0.1;...% Thalumus, light brown
                 'Dropout',256,0.1,0.1,0.1}; % Signal dropout regions
             
for ii = 1:length(module_id_rgb)
    ind = module_id_rgb{ii,2};
    cmap_base(ind,:) = cell2mat(module_id_rgb(ii,3:end));
end

% Read in data (vector of values corresponding to each ROI in the label/surface files)
% Column format
data = dlmread('test_data_integer.txt');

% Color range (min val, max val)
colorrange = [1,256];

% Surface gifti file
fs_dir = 'fsaverage_BrainnetomeAtlas/164k/' ;

% Surface gifti file
surface = [fs_dir 'fsaverage.L.very_inflated.164k.surf.gii'];

% Label gifti file
label = [fs_dir 'fsaverage.L.BN_Atlas.164k.label.gii'];

% Output directory
out_dir = pwd;

% Output file name
out_names = {'Surf_Integer_lat_l_270_0';...
             'Surf_Integer_med_l_90_0';...
             'Surf_Integer_sup_l_0_90'};

% Viewing angle
view_angs = [270,0;...
            90,0;...
            0,90];

% ----------- RUN FUNCTION --------------
for ii = 1:numel(out_names)
    out_name = out_names{ii};
    view_ang = view_angs(ii,:);
    plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);
end