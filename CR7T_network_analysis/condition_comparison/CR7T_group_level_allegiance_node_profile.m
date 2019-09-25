%% CR7T - Test node level allegiance profile between conditions (at group level)
% This script assesses the similarity of node level allegiance profiles
% between the symbolic and nonsymbolic conditions. Similarity is defined as
% the Pearson correlation between the node's allegiance vector for each
% condition, i.e. the pattern of allegiance from the node to all other
% nodes. The Fisher Z value is then tested against the non-parametric null
% distribution, to assess whether the degree of (dis)similarity is significant.
% Ultimately, we are looking for nodes which show a significantly different
% pattern of allegiance between formats..

purge
% Load allegiance matrices
load('workspace_consensus_clustering_gamma_2pt45','Pa_group');
% Load allegiance matrices from permutations
load('workspace_group_consensus_pair_permutations_50000iters_03_27_05_13_22.mat','Pa_group_perm','num_perms');

% Load module assignments and convert to 246 vector
load('workspace_Qc_significance_symbolic.mat', 'C');
C_sym = CR7T_convert_vec_BN202_to_BN246(C);
load('workspace_Qc_significance_nonsymbolic.mat', 'C');
C_nonym = CR7T_convert_vec_BN202_to_BN246(C);

%% Calculate similarity of the node allegiance profile between conditions
Pa_s = Pa_group.BN202_Compare_Symbolic;
Pa_n = Pa_group.BN202_Compare_Nonsymbolic;
num_nodes = size(Pa_s,1);
for ii = 1:num_nodes
    r = corr(Pa_s(:,ii),Pa_n(:,ii));
    r_true(ii) = r;
end
% Transform to Fisher Z
z_true = atanh(r_true);

%% Calculate null distributions of similarity values
for ii = 1:num_perms
    disp(num2str(ii));
    Pa1 = Pa_group_perm(:,:,ii,1);
    Pa2 = Pa_group_perm(:,:,ii,2);
    for jj = 1:num_nodes
        r = corr(Pa1(:,jj),Pa2(:,jj));
        r_perm(jj,ii) = r;
    end
end
% Transform to Fisher Z
z_perm = atanh(r_perm);

%% Get significance
for ii = 1:num_nodes
    zscore_true(ii) = (z_true(ii) - mean(z_perm(ii,:))) / std(z_perm(ii,:));
    pval_true(ii) = invprctile(z_perm(ii,:),z_true(ii));
end

%% Transpose vectors into combined results table
results(:,1) =  CR7T_convert_vec_BN202_to_BN246(r_true);
results(:,2) =  CR7T_convert_vec_BN202_to_BN246(z_true);
results(:,3) =  CR7T_convert_vec_BN202_to_BN246(zscore_true);
results(:,4) =  CR7T_convert_vec_BN202_to_BN246(pval_true);

%% Save workspace
clear Pa_group_perm
save('node_profile_statistics.mat');

%% Write edge matrix for the top nodes (Pa difference)
% Get difference matrix
Pa_diff = Pa_s - Pa_n;
% Get sorted order of nodes based on
nodes = 1:length(Pa_diff);
[~,sortind] = sort(zscore_true);
nodes_sorted = nodes(sortind);

% Load node info
load('workspace_labels.mat','T');
labels = T.Index_new;

% Prepare matrix for BN202 node coordinate files
out_dir = 'Z:/CR7T_Connectivity/matfiles/CR7T_network_analysis/condition_comparison/BrainNetViewer_Plotting/';
BN246_nodefile = load('BN246_MNIcoord4BrainNetViewer.node');
BN202_nodefile = BN246_nodefile;
BN202_nodefile(labels == 0,:) = [];
BN202_nodefile(:,6) = nodes;
dlmwrite('BN202_MNIcoord4BrainNetViewer.node',BN202_nodefile,'delimiter',' ');

% Make edge file for each of the top nodes
num_sig = sum(pval_true < 5); 
for ii = 1:num_sig
    % Node index in BN202
    node = nodes_sorted(ii);
    nvec = Pa_diff(node,:);
    % Set up matrix for only that node
    nmat = zeros(size(Pa_diff,1),size(Pa_diff,2));
    nmat(node,:) = nvec;
    nmat(:,node) = nvec';
    % Get node label
    nlabel = T.Region(labels == node);
    
    % Create positive, negative, and combined thresholded matrix
    thresh = 0.12; % With 33 subjects, this indicates different in allegiance of 4 or more subjects
    nmat_pos = nmat;
    nmat_pos(nmat_pos < thresh) = 0;
    nmat_neg = nmat;
    nmat_neg(nmat_neg > -thresh) = 0;
    nmat_combine = nmat_pos + nmat_neg;
    
    % Create node file with region of interest as module 1, 
    nodefile = BN202_nodefile;
    % Get indices of surviving positive and negative weights
    pos = nmat_pos(node,:) > 0;
    neg = nmat_neg(node,:) < 0;
    % Set "module" label of nodes, for coloring
    nodefile(node,4) = 1;
    nodefile(pos,4) = 2;
    nodefile(neg,4) = 3;
    % Set size of nodes (proportional to the weight of the difference,
    % after some scaling)
    nodefile(:,5) = 3*(nmat_pos(node,:).^2 + abs(nmat_neg(node,:)).^2);
%     nodefile(:,5) = nmat_pos(node,:).*0.75 + abs(nmat_neg(node,:)).*0.75;
%     nodefile(:,5) = nmat_pos(node,:) + abs(nmat_neg(node,:));
    nodefile(node,5) = 0.5;
    
    % Write files
    prefix = [out_dir 'Pa_Diff_Sym-Nonsym_Gamma_2.45.' nlabel{1}];
    dlmwrite([out_dir 'BN202_' nlabel{1} '_' num2str(thresh) '.node'],nodefile,'delimiter',' ');
    dlmwrite([prefix '_pos_' num2str(thresh) '.edge'],nmat_pos,'delimiter',' ');
    dlmwrite([prefix '_neg_' num2str(thresh) '.edge'],nmat_neg,'delimiter',' ');
    dlmwrite([prefix '_combine_' num2str(thresh) '.edge'],nmat_combine,'delimiter',' ');
    dlmwrite([prefix '.edge'],nmat,'delimiter',' ');
end

%% Vertical bar chart of top values
figure('Position',[100,100,400,1000]);
zscore_true_sig = zscore_true(sortind);
zscore_true_sig = flip(zscore_true_sig(1:num_sig));
b = barh(abs(zscore_true_sig));
b = barh(zscore_true_sig);
b.BarWidth = .7;
grid on
ax = gca;
ax.XLabel.String = 'Z-score';
for ii = 1:num_sig
    node = nodes_sorted(ii);
    ax.YTickLabel(ii) =  T.Region(T.Index_new == node);
end
ax.YTickLabel = flip(ax.YTickLabel);
ax.TickLabelInterpreter = 'none';
ax.FontSize = 15;
ax.XTick = [0 0.5 1 1.5 2 2.5 3];
ax.XTickLabel = {'0','-0.5','-1','-1.5','-2','-2.5','-3'};
ax.Title.String = 'Allegiance Dissimilarity (p<0.05)';


figure('Position',[100,100,600,1000]);
zscore_true_sig = zscore_true(sortind);
zscore_true_sig = flip(zscore_true_sig(1:num_sig));
b = barh(zscore_true_sig);
b.BarWidth = .7;
grid on
ax = gca;
ax.XLabel.String = 'Z-score';
for ii = 1:num_sig
    node = nodes_sorted(ii);
    ax.YTickLabel(ii) =  T.Region(T.Index_new == node);
end
ax.YTickLabel = flip(ax.YTickLabel);
ax.TickLabelInterpreter = 'none';
ax.FontSize = 15;
%ax.XTick = [-3 -2.5 -2 -1.5 -1 -.5 0];
%ax.XTickLabel = {'0','-0.5','-1','-1.5','-2','-2.5','-3'}; ax.XTickLabel = flip(ax.XTickLabel);
ax.Title.String = 'Allegiance Dissimilarity (p<0.05)';
ax.XLim = [-3 -1.5];

%% Write niftis with node level similarity values
% % Set output suffix
% out_dir = '/mnt/CR7T_Connectivity/matfiles/CR7T_network_analysis/condition_comparison/';
% % Set atlas base
% atlas_nii = '/mnt/CR7T_Connectivity/atlas/BN_Atlas_246_1mm.nii.gz';
% % Set out names 
% out_name = {'Allegiance_profile_similarity_Sym_vs_Nonsym_FisherZ';...
%             'Allegiance_profile_similarity_Sym_vs_Nonsym_Zscore_perm';...
%             'Allegiance_profile_similarity_Sym_vs_Nonsym_Zscore_perm_p.05';...
%             'Allegiance_profile_similarity_Sym_vs_Nonsym_Zscore_perm_p.01';...
%             'Allegiance_profile_similarity_Sym_vs_Nonsym_Zscore_perm_p.005'};
% % Set data for each nifti
% c2(:,1)= CR7T_convert_vec_BN202_to_BN246(z_true);
% c2(:,2)= CR7T_convert_vec_BN202_to_BN246(zscore_true);
% c = zscore_true; c(pval_true > 5) = 0; % P < 0.05
% c2(:,3) = CR7T_convert_vec_BN202_to_BN246(c);
% c = zscore_true; c(pval_true > 1) = 0; % P < 0.01
% c2(:,4) = CR7T_convert_vec_BN202_to_BN246(c);
% c = zscore_true; c(pval_true > 0.5) = 0; % P < 0.005
% c2(:,5) = CR7T_convert_vec_BN202_to_BN246(c);
%             
% for ii = 1:numel(out_name)
%     % Nifti with all modules
%     txt = [out_dir out_name{ii} '.txt'];
%     nii = [out_dir out_name{ii} '.nii.gz'];
%     dlmwrite(txt,c2(:,ii),'delimiter',' ');
%     unix(['rm -f ' nii]);
%     unix(['3dUndump -datum float -ROImask ' atlas_nii ' -prefix ' nii ' ' txt]);
% end