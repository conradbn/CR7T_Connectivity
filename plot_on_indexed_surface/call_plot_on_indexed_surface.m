%% CALL SURFACE PLOTTING FUNCTION
% Benjamin Neal Conrad, Vanderbilt University, November 2018

clear all; close all

%% Read in data (vector of values corresponding to each ROI in the label/surface files)
data = dlmread('test_data.txt');
data = data(1,:);

%% Base colormap
% % RdBu colormap (requires cbrewer package)
% RdBu = cbrewer('div', 'RdBu', 256);
% c = flip(RdBu); 
c = jet(256);
cmap_base = c;

%% Color range (min val, max val)
colorrange = [-1.5,1.5];

%% Surface gifti file
fs_dir = 'fsaverage_BrainnetomeAtlas\164k\' ;
surface = [fs_dir 'fsaverage.R.very_inflated.164k.surf.gii'];

%% Label gifti file
label = [fs_dir 'fsaverage.R.BN_Atlas.164k.label.gii'];

%% Output directory
out_dir = pwd;

%% Output file name
out_name = 'test.png';

%% Viewing angle
view_ang = [90,0];

%% ----------- RUN FUNCTION --------------
plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang);

