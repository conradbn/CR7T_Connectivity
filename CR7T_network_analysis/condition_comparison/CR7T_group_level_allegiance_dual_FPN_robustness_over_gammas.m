%% Check the robustness of dual FPN finding over gammas
purge;
% Load partitions
load('workspace_consensus_clustering_gamma_2pt45.mat','C_group'); 
cnames = fieldnames(C_group);
% Relabel partitions
C_relabel = C_group;
C_relabel.(cnames{1}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{1}));
C_relabel.(cnames{2}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{2}));
C_sym = C_relabel.BN202_Compare_Symbolic;
C_nonsym = C_relabel.BN202_Compare_Nonsymbolic;

load('workspace_modularity_gamma_sweep_33subs_Symbolic_1000xLouvain_1xGamma_2pt45.mat','C_group_iter','gammas');
C_sym_sweep = C_group_iter;
load('workspace_modularity_gamma_sweep_33subs_Nonsymbolic_1000xLouvain_1xGamma_2pt45.mat','C_group_iter');
C_nonsym_sweep = C_group_iter;

%% Relabel each partition
for ii = 1:size(C_sym_sweep,2)
    C_sym_sweep(:,ii) = pair_labeling(C_sym,C_sym_sweep(:,ii));
    C_nonsym_sweep(:,ii) = pair_labeling(C_nonsym,C_nonsym_sweep(:,ii));
    %C_nonsym_sweep(:,ii) = pair_labeling(C_sym,C_nonsym_sweep(:,ii));
end

%% Count number of "VAN" and "DAN" regions in each partition
for ii = 1:size(C_sym_sweep,2)
    count_sym(ii,1) = sum(C_sym_sweep(:,ii) == 6);
    count_sym(ii,2) = sum(C_sym_sweep(:,ii) == 24);
    count_nonsym(ii,1) = sum(C_nonsym_sweep(:,ii) == 6);
    count_nonsym(ii,2) = sum(C_nonsym_sweep(:,ii) == 24);
end

%% Get all relevant nodes/labels
ind = C_sym == 6 ;
ind = or(ind,C_sym == 24);
ind = or(ind,C_nonsym == 6);
nodes_sym = C_sym_sweep(ind,:);
nodes_nonsym = C_nonsym_sweep(ind,:);
nodes_sym = sort(nodes_sym);
nodes_nonsym = sort(nodes_nonsym);

% Get colormap (Sym has more unique)
uniq = unique(nodes_sym);
uniqns = unique(nodes_nonsym);
nuniq = numel(uniq);
cmap = distinguishable_colors(nuniq);
% Set the values for DAN/VAN/FPN
cmap(uniq == 6, :) = [1,1,0]; %fpn/van
cmap(uniq == 24, :) = [0,1,0]; %dan
cmap(uniq == 2, :) = [1,0,0]; %dmn
cmap(uniq == 21, :) = [0,0,1]; %visual
cmap(uniq == 35, :) = [0.8,.33,0]; %dorsal caudate
% Translate colors into new values
nodes_sym2 = nodes_sym;
nodes_nonsym2 = nodes_nonsym;
for ii = 1:numel(uniq)
    nodes_sym2(nodes_sym == uniq(ii)) = ii;
end
for ii = 1:numel(uniqns)
    nodes_nonsym2(nodes_nonsym == uniqns(ii)) = ii;
end
  
% NONSYMBOLIC
figure('Position',[100,100,900,600]);
imagesc(nodes_nonsym2);
% nuniq = numel(unique(nodes_nonsym));
% cmap = distinguishable_colors(nuniq);
colormap(cmap);
ax = gca;
%ax.CLim = [0 nuniq];
ax.Title.String = 'Nonsymbolic - FPN/DAN/VAN Node Assignments';
ax.Title.FontSize = 17;
% X Ticks
ax.XTick = 0.5:numel(gammas);
gcells = num2cell(gammas(1:end));
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 260;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
% Axis labels
ax.XLabel.String = 'Gamma Value';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 15;
ax.YLabel.String = 'Nodes (arbitrary order)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 15;
ax.YTick = [];
grid on

% SYMBOLIC
figure('Position',[100,100,900,600]);
imagesc(nodes_sym2);
colormap(cmap);
ax = gca;
%ax.CLim = [0 nuniq];
ax.Title.String = 'Symbolic - FPN/DAN/VAN Node Assignments';
ax.Title.FontSize = 17;
%grid on;
% X Ticks
ax.XTick = 0.5:numel(gammas);
gcells = num2cell(gammas(1:end));
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 260;
ax.XAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
% Axis labels
ax.XLabel.String = 'Gamma Value';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 15;
ax.YLabel.String = 'Nodes (arbitrary order)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 15;
ax.YTick = [];
grid on


%% Plot
figure;
plot(count_nonsym)
hold on; 
plot(count_sym)
grid on
ax = gca;
ax.FontSize = 17;
%ax.YLim = [0.1 0.5];
%ax.GridAlpha = .3;
% Title
ax.Title.String = 'Variability of Flexibility Across Nodes - Gamma Sweep';
ax.Title.Interpreter = 'none';
% X Ticks
ax.XTick = 1:2:numel(gammas);
gcells = num2cell(gammas(1:2:end));
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 270;
ax.XAxis.FontSize = 12;
 ax.XAxis.FontWeight = 'bold';
% Axis labels
xL = xlabel('Gamma Value');
xL.FontWeight = 'bold';
xL.FontSize = 12;
ax.YLabel.String = 'Number of VAN/DAN nodes over gammas';



