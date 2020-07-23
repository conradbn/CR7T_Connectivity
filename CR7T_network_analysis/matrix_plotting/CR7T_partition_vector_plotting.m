%% Consensus partition vectors
purge

load('Z:\CR7T_Connectivity\CR7T_Connectivity-Git\CR7T_data_files\workspace_consensus_clustering_gamma_2pt45.mat','M_all','C_group','C_all','Pa_group'); 
%load('Z:\CR7T_Connectivity\matfiles\CR7T_WORKSPACES\workspace_consensus_clustering_gamma_2pt45.mat','M_all','C_group','C_all','Pa_group'); 
cnames = fieldnames(C_group);
C_relabel = C_group;
C_relabel.(cnames{1}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{1}));
C_relabel.(cnames{2}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{2}));

%% Reording
% Reorder using all trials conditon

% All
Pa = Pa_group.BN202_Compare_All;
Cm = C_relabel.BN202_Compare_All;
% Perform reordering function
[On,~] = reorder_mod(Pa,Cm);

% Symbolic
Pa = Pa_group.BN202_Compare_Symbolic;
Cm = C_relabel.BN202_Compare_Symbolic;
% Perform reordering function
[On_s,~] = reorder_mod(Pa,Cm);

% Nonsymbolic
Pa = Pa_group.BN202_Compare_Nonsymbolic;
Cm = C_relabel.BN202_Compare_Nonsymbolic;
% Perform reordering function
[On_n,~] = reorder_mod(Pa,Cm);

% % Get module box lines
% c = Cm(On);
% nc = max(c);
% X = [];
% Y = [];
% val_order = unique(c,'stable');
% 
% for i = 1:nc
%     ind = find(c == val_order(i));
%     if ~isempty(ind)
%         mn = min(ind) - 0.5;
%         mx = max(ind) + 0.5;
%         x = [mn mn mx mx mn NaN];
%         y = [mn mx mx mn mn NaN];
%         X = [X, x]; 
%         Y = [Y, y];
%     end
% end

%% Final group-level module colormap
module_id_rgb = {'DMN',2,1,0,0;... % DMN, red
                 'FPN',6,1,1,0;... % FPN, yellow
                 'SMN',10,0,1,1;... % Sensorimotor, cyan
                 'Sal',12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
                 'STS',16,1,0.7,1;... % STS, pink
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
ind = find(all(cmap_base==0,2));
cmap_base(ind,:) = 0.7;

%% Set subject partitions
s = 33;
ss = [1,16,33];% 12 26
MsMnCsCn = [squeeze(M_all.BN202_Compare_Symbolic(:,1:3,1,s)),...
            squeeze(M_all.BN202_Compare_Nonsymbolic(:,1:3,1,s)),...
            squeeze(C_all.BN202_Compare_Symbolic(:,1,ss)),...
            squeeze(C_all.BN202_Compare_Nonsymbolic(:,1,ss))];

%% Plot the subject partitions
%base_s = pair_labeling(C_relabel.BN202_Compare_Nonsymbolic,C_relabel.BN202_Compare_Symbolic);
base_s = C_relabel.BN202_Compare_Symbolic;
base_n = C_relabel.BN202_Compare_Nonsymbolic;
colors = colorcube(numel(unique(MsMnCsCn(:,8))));
for m = 1:size(MsMnCsCn,2)
    if m<=3 || (m>6 && m<10)
        v = pair_labeling(base_s,MsMnCsCn(:,m));
        v = v(On);
    else
        v = pair_labeling(base_n,MsMnCsCn(:,m));
        v = v(On);
    end

    figure('Position',[3000 100 100 1000]); hold on;
    axis([0 3 0 202])
    for ii = numel(v):-1:1
        node = v(ii);
        %color = colors(v(ii),:);
        color = cmap_base(v(ii),:);
        pos = 202-ii;
        rectangle('Position',[1 pos 1 1],'FaceColor',color);
    end
    axis off
    %print(['MsMnCsCn' num2str(m) '.png'],'-dpng')
    export_fig(['MsMnCsCn' num2str(m) '.png'], '-dpng','-transparent');
end

%% Plot the final group partitions
v = C_relabel.BN202_Compare_Nonsymbolic;
v = v(On_n);
v(v==16|v==20|v==24) = 1;
figure('Position',[3000 100 100 1000]); hold on;
axis([0 3 0 202])
for ii = numel(v):-1:1
    node = v(ii);
    %color = colors(v(ii),:);
    color = cmap_base(v(ii),:);
    pos = 202-ii;
    rectangle('Position',[1 pos 1 1],'FaceColor',color);
end
axis off
export_fig('Nonsymbolic_partition_reordered.png', '-dpng','-transparent');



v = C_relabel.BN202_Compare_Symbolic;
v = v(On_s);
v(v==12|v==27|v==35) = 1;
figure('Position',[3000 100 100 1000]); hold on;
axis([0 3 0 202])
for ii = numel(v):-1:1
    node = v(ii);
    %color = colors(v(ii),:);
    color = cmap_base(v(ii),:);
    pos = 202-ii;
    rectangle('Position',[1 pos 1 1],'FaceColor',color);
end
axis off
export_fig('Symbolic_partition_reordered.png', '-dpng','-transparent');








