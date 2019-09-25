%% Plot group level allegiance matrices, reordered, and difference in conditions
purge
% Load data 
load('workspace_consensus_clustering_gamma_2pt45.mat','Pa_group','C_group');
cnames = fieldnames(C_group);
% Relabel partitions
C_relabel = C_group;
C_relabel.(cnames{1}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{1}));
C_relabel.(cnames{2}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{2}));

%% SYMBOLIC
% Set matrices of interest
Pa = Pa_group.BN202_Compare_Symbolic;
Cm = C_relabel.BN202_Compare_Symbolic;

% Reording
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

% Plotting
% Get colormap
RdBu = cbrewer('div', 'RdBu', 256);
RdBu = flip(RdBu); 
% Plot matrix
fig = figure('Position',[100,100,900,900]);
imagesc(Pa(On,On)); colorbar; colormap(RdBu);
ax = gca; 
ax.Title.String = 'Symbolic Allegiance Matrix - Symbolic Modules';
ax.FontSize = 18;
ax.CLim = [0 0.7];
ax.XTick = [];
ax.YTick = [];
hold on;
plot(X,Y,'k','Linewidth',2);
axis square; hold off;
out = strrep(ax.Title.String,' ','_');
print(fig,out,'-dpng');

% Plot the other Pa matrix
Pa2 = Pa_group.BN202_Compare_Nonsymbolic;
fig = figure('Position',[100,100,900,900]);
imagesc(Pa2(On,On)); colorbar; colormap(RdBu);
ax = gca; 
ax.Title.String = 'Nonsymbolic Allegiance Matrix - Symbolic Modules';
ax.FontSize = 18;
ax.CLim = [0 0.7];
ax.XTick = [];
ax.YTick = [];
hold on;
plot(X,Y,'k','Linewidth',2);
axis square; hold off;
out = strrep(ax.Title.String,' ','_');
print(fig,out,'-dpng');

% Plot the Difference matrix
Pa_diff = Pa - Pa2;
fig = figure('Position',[100,100,900,900]);
imagesc(Pa_diff(On,On)); colorbar; colormap(RdBu);
ax = gca; 
ax.Title.String = 'Allegiance Matrix Difference - Symbolic Modules';
ax.FontSize = 18;
ax.CLim = [-0.25 0.25];
ax.XTick = [];
ax.YTick = [];
hold on;
plot(X,Y,'k','Linewidth',2);
axis square; hold off;
out = strrep(ax.Title.String,' ','_');
print(fig,out,'-dpng');

%% NONSYMBOLIC
% Set matrices of interest
Pa = Pa_group.BN202_Compare_Nonsymbolic;
Cm = C_relabel.BN202_Compare_Nonsymbolic;

% Reording
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

% Plotting
% Get colormap
RdBu = cbrewer('div', 'RdBu', 256);
RdBu = flip(RdBu); 
% Plot matrix
fig = figure('Position',[100,100,900,900]);
imagesc(Pa(On,On)); colorbar; colormap(RdBu);
ax = gca; 
ax.Title.String = 'Nonsymbolic Allegiance Matrix - Nonsymbolic Modules';
ax.FontSize = 18;
ax.CLim = [0 0.7];
ax.XTick = [];
ax.YTick = [];
hold on;
plot(X,Y,'k','Linewidth',2);
axis square; hold off;
out = strrep(ax.Title.String,' ','_');
print(fig,out,'-dpng');

% Plot the other Pa matrix
Pa2 = Pa_group.BN202_Compare_Symbolic;
fig = figure('Position',[100,100,900,900]);
imagesc(Pa2(On,On)); colorbar; colormap(RdBu);
ax = gca; 
ax.Title.String = 'Symbolic Allegiance Matrix - Nonsymbolic Modules';
ax.FontSize = 18;
ax.CLim = [0 0.7];
ax.XTick = [];
ax.YTick = [];
hold on;
plot(X,Y,'k','Linewidth',2);
axis square; hold off;
out = strrep(ax.Title.String,' ','_');
print(fig,out,'-dpng');

% Plot the Difference matrix
Pa_diff = Pa - Pa2;
fig = figure('Position',[100,100,900,900]);
imagesc(Pa_diff(On,On)); colorbar; colormap(RdBu);
ax = gca; 
ax.Title.String = 'Allegiance Matrix Difference - Nonsymbolic Modules';
ax.FontSize = 18;
ax.CLim = [-0.25 0.25];
ax.XTick = [];
ax.YTick = [];
hold on;
plot(X,Y,'k','Linewidth',2);
axis square; hold off;
out = strrep(ax.Title.String,' ','_');
print(fig,out,'-dpng');
