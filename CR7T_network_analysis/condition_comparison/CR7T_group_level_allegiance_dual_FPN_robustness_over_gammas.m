%% Check the robustness of dual FPN finding over gammas
% Edited 1/24/2020 - larger gamma range and edited community names
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

%load('workspace_modularity_gamma_sweep_33subs_Symbolic_1000xLouvain_1xGamma_2pt45.mat','C_group_iter','gammas');
%load('workspace_modularity_gamma_sweep_33subs_Symbolic_1000xLouvain_1xGamma_0pt5_to_3pt5_savePa_Call.mat','C_group_iter');
load('workspace_modularity_gamma_sweep_33subs_Symbolic_1000xLouvain_1xGamma_0_to_5.mat','C_group_iter','gammas');
C_sym_sweep = C_group_iter(:,2:end);
gammas = gammas(2:end);
%load('workspace_modularity_gamma_sweep_33subs_Nonsymbolic_1000xLouvain_1xGamma_2pt45.mat','C_group_iter');
%load('workspace_modularity_gamma_sweep_33subs_Nonsymbolic_1000xLouvain_1xGamma_0pt5_to_3pt5_savePa_Call.mat','C_group_iter','gammas');
load('workspace_modularity_gamma_sweep_33subs_Nonsymbolic_1000xLouvain_1xGamma_0_to_5.mat','C_group_iter');
C_nonsym_sweep = C_group_iter(:,2:end);

%% Relabel each partition
for ii = 1:size(C_sym_sweep,2)
    C_sym_sweep(:,ii) = pair_labeling(C_sym,C_sym_sweep(:,ii));
    C_nonsym_sweep(:,ii) = pair_labeling(C_nonsym,C_nonsym_sweep(:,ii));
    %C_sym_sweep(:,ii) = pair_labeling(C_nonsym,C_sym_sweep(:,ii));
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
ind = C_sym == 6;
ind = or(ind,C_sym == 24);
ind = or(ind,C_nonsym == 6);
nodes_sym = C_sym_sweep(ind,:);
nodes_nonsym = C_nonsym_sweep(ind,:);

%% Relabel again to ensure maximal consistency in colors across both plots
% %vals_to_check = [2,6,10,16,20,21,24,35,36];
% s = nodes_sym(:,gammas == 2.45);
% n = nodes_sym(:,gammas == 2.45);
% for g = 1:size(nodes_sym,2)
%     nodes_sym(:,g) = pair_labeling(s,nodes_sym(:,g));
%     nodes_nonsym(:,g) = pair_labeling(n,nodes_nonsym(:,g));
% end

%% Plotting
% Sort nodes 
nodes_sym = sort(nodes_sym);
nodes_nonsym = sort(nodes_nonsym);
% Get colormap (Sym has more unique)
uniq = unique(nodes_sym);
uniqns = unique(nodes_nonsym);
nuniq = numel(uniq);
cmap = distinguishable_colors(nuniq);
% Set the values for DAN/FPN
cmap(uniq == 6, :) = [1,1,0]; %fpn
cmap(uniq == 24, :) = [0,1,0]; %dan
cmap(uniq == 2, :) = [1,0,0]; %dmn
cmap(uniq == 21, :) = [0,0,1]; %visual
cmap(uniq == 10, :) = [0,1,1]; % SMN
cmap(uniq == 16, :) = [1,0.7,1];% Aud
cmap(uniq == 20, :) = [1,0.5,0];% CdV/NAc
cmap(uniq == 35, :) = [.8,.33,0]; % CdD
cmap(uniq == 36, :) = [0.5,0.25,0.1]; % Thal

% Translate colors into new values
nodes_sym2 = nodes_sym;
nodes_nonsym2 = nodes_nonsym;
for ii = 1:numel(uniq)
    nodes_sym2(nodes_sym == uniq(ii)) = ii;
end
for ii = 1:numel(uniqns)
    nodes_nonsym2(nodes_nonsym == uniqns(ii)) = ii;
end

% % Remove black color
% cmap(sum(cmap,2)<0.18,:) = [0.4,0.6,0.4];

% NONSYMBOLIC
figure('Position',[100,100,1400,500]);
imagesc(nodes_nonsym2);
% nuniq = numel(unique(nodes_nonsym));
% cmap = distinguishable_colors(nuniq);
colormap(cmap);
ax = gca;
%ax.CLim = [0 nuniq];
ax.Title.String = 'Nonsymbolic - FPN/DAN Region Assignments';
ax.Title.FontSize = 25;
% X Ticks
ax.XTick = 0.5:numel(gammas);
gcells = num2cell(gammas(1:end));
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 270;
ax.XAxis.FontSize = 10;
ax.XAxis.FontWeight = 'bold';
% Axis labels
ax.XLabel.String = 'Gamma Value (0.05 - 5)';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 20;
ax.YLabel.String = 'Regions (arbitrary order)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 20;
ax.YTick = [];
grid on
%print('FPN-DAN_region_assignments_gamma_sweep_Nonsymbolic.png','-dpng');
export_fig('FPN-DAN_region_assignments_gamma_sweep_Nonsymbolic','-png','-m2','-transparent');
% SYMBOLIC
figure('Position',[100,100,1400,500]);
imagesc(nodes_sym2);
colormap(cmap);
ax = gca;
%ax.CLim = [0 nuniq];
ax.Title.String = 'Symbolic - FPN/DAN Region Assignments';
ax.Title.FontSize = 25;
%grid on;
% X Ticks
ax.XTick = 0.5:numel(gammas);
gcells = num2cell(gammas(1:end));
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 270;
ax.XAxis.FontSize = 10;
ax.XAxis.FontWeight = 'bold';
% Axis labels
ax.XLabel.String = 'Gamma Value (0.05 - 5)';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 20;
ax.YLabel.String = 'Regions (arbitrary order)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 20;
ax.YTick = [];
grid on
%print('FPN-DAN_region_assignments_gamma_sweep_Symbolic.png','-dpng');
export_fig('FPN-DAN_region_assignments_gamma_sweep_Symbolic','-png','-m2','-transparent');


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
ax.YLabel.String = 'Number of FPN/DAN nodes over gammas';



