%% Reorder the connectivity z-score matrices
purge
% Load data 
load('workspace_consensus_clustering_gamma_2pt45.mat');
load('Z:/CR7T_Connectivity/matfiles/CR7T_WORKSPACES/workspace_zscore_mats_33subs.mat');

cnames = fieldnames(C_group);
% Reorder using all trials conditon
% Set matrices of interest
Pa = Pa_group.(cnames{2});
Cm = C_group.(cnames{2});

%% Reording
% Perform reordering function
[On,Wr] = reorder_mod(Pa,Cm);

% Get module box lines
c = Cm(On);
nc = max(c);
X = [];
Y = [];
val_order = unique(c,'stable');

for i = 1:nc
    ind = find(c == val_order(i));
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        mx = max(ind) + 0.5;
        x = [mn mn mx mx mn NaN];
        y = [mn mx mx mn mn NaN];
        X = [X, x]; 
        Y = [Y, y];
    end
end

%% Plotting
% Get colormap
RdBu = cbrewer('div', 'RdBu', 256);
RdBu = flip(RdBu); 
% Plot matrix
figure('Position',[100,100,900,900]);
imagesc(Pa(On,On)); colorbar; colormap(RdBu); %ax = gca; colormap(ax,cmap_081517); colorbar; %caxis([0 .5]);
hold on;
plot(X,Y,'k','Linewidth',3);
title([cnames{1} ' Allegiance Matrix - Reordered by Module'],'interpreter','none'); 
axis square; hold off;

% Do the conditions separate
for cond = 2:numel(cnames)
    % Set matrices of interest
    Pa = Pa_group.(cnames{cond});
    Cm = C_group.(cnames{cond});
    % Plot
    figure('Position',[100,100,900,900]);
    imagesc(Pa(On,On)); colorbar; colormap(RdBu); %ax = gca; colormap(ax,cmap_081517); colorbar; %caxis([0 .5]);
    hold on;
    plot(X,Y,'k','Linewidth',3);
    title([cnames{cond} ' Allegiance Matrix - Reordered by Module'],'interpreter','none'); 
    axis square; hold off;
end

%% Plot all subject allegiance matrices
for cond = 1:numel(cnames)
    % Load condition mats, remove inf/nan vals
    Pa_subj = squeeze(Pa_all.(cnames{cond}));
    % Set figure window size for better looking display
    figure('Position', [100, 100, 1700, 1000]);
    string = strsplit(cnames{cond},'_');
    % Plot matrices
    for jj = 1:size(Pa_subj,3)
        %subaxis(5,8,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
        subaxis(5,8,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
        imagesc(Pa_subj(On,On,jj)); colormap(RdBu); 
        caxis([0 1])
        axis square;
        axis off; 
        subj_name = subj_dirs(jj).name;
        title([string{3} ' Allegiance CR-' subj_name(4:6)],'interpreter','none');
        set(gca,'FontSize',8);
    end
end

%% Plot difference Matrix
% Plot matrix
Pa_diff = Pa_group.BN202_Compare_Symbolic - Pa_group.BN202_Compare_Nonsymbolic;
figure('Position',[100,100,900,900]);
imagesc(Pa_diff(On,On)); colorbar; colormap(RdBu); %ax = gca; colormap(ax,cmap_081517); colorbar; %caxis([0 .5]);
hold on;
plot(X,Y,'k','Linewidth',3);
title([cnames{1} ' Allegiance Difference Matrix - Reordered by Module'],'interpreter','none'); 
axis square; hold off;

%% Connectivity matrices
%close all
% load('zscore_mats_102018_final.mat')
% Plot all subject connectivity matrices
for cond = 1:numel(cnames)
    % Load condition mats, remove inf/nan vals
    zs_subj = squeeze(zscore_mats.(cnames{cond}));
    % Set figure window size for better looking display
    figure('Position', [100, 100, 1700, 1000]);
    string = strsplit(cnames{cond},'_');
    % Plot matrices
    for jj = 1:size(zs_subj,3)
        %subaxis(5,8,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
        subaxis(5,8,jj,'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05); 
        imagesc(zs_subj(On,On,jj)); colormap(RdBu); 
        caxis([-8 8])
        axis square;
        axis off; 
        subj_name = subj_dirs(jj).name;
        title([string{3} ' Connectivity CR-' subj_name(4:6)],'interpreter','none');
        set(gca,'FontSize',7);
    end
end