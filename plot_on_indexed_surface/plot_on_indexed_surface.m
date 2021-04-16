function out = plot_on_indexed_surface(data,cmap_base,colorrange,surface,label,out_dir,out_name,view_ang)
%% Plot values to ROIs on surface of gifti images
% Benjamin Neal Conrad, Vanderbilt University, November 2018

%% Map color range to colormap
crange_map = linspace(colorrange(1),colorrange(2),length(cmap_base));
%% Load Surface file
g = gifti(surface);
%% Load Label file
gg = gifti(label);
ggl = gg.labels;
num_rois = numel(ggl.name);

%% For each vector of data
% Setup up new color map with number of rois
cmap = cmap_base(1:num_rois,:);
% Transpose data
inds = data;
% Insert zero for label 1
inds = [0;inds]';
% Loop through each ROI
for jj = 1:size(cmap,1)
    % Get value for ROI
    ind = inds(jj);
    % Set values of 0 to gray
    if ind == 0
        cmap(jj,:) = [.6,.6,.6]; % Hardcoded to "gray"
    else
        % Set color for ROI based on input colormap
        [~,idx]=min(abs(crange_map-ind));
        minVal=crange_map(idx);
        cmap(jj,:) = cmap_base(idx,:);
    end
end

%% Plot
fig = figure;
set(gcf,'visible','off');
h = plot(g,gg);
view(view_ang)
lighting flat % gouraud
lightangle(90,20)
material dull
colormap(cmap);
% Prepare fig for printing and save as png
out = [out_dir '/' out_name];
export_fig(out,'-transparent','-r300');
