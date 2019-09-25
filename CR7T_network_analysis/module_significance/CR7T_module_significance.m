%% CR7T - Statistical Significance of Modules
%% Set variables
purge
% Set condition of interest
cond_modules = 'BN202_Compare_Symbolic';
cond_data = 'BN202_Compare_Symbolic';
% Set number of permutations
nperms = 10000;
% Community louvain parameters
gamma = 2.45;
tau = 0.5;
B = 'negative_asym';

%% Load/setup data
% Load consensus partitions
load('workspace_consensus_clustering_gamma_2pt45','zscore_mats','C_group');
cnames = fieldnames(C_group);
% Relabel partitions
C_relabel = C_group;
C_relabel.(cnames{1}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{1}));
C_relabel.(cnames{2}) = pair_labeling(C_relabel.(cnames{3}),C_relabel.(cnames{2}));

% Get consensus partition of interest
C = C_relabel.(cond_modules);
% Get the unique module indices
modules_all = unique(C);

% Remove modules if singleton
ind = 1;
for ii = 1:numel(modules_all)
    module = modules_all(ii);
    if nnz(C==module) ~= 1
        modules_true(ind) = module;
        ind = ind + 1;
    end
end

%% Calculate modularity Q for each module
% Connectivity matrices
mats = zscore_mats.(cond_data);
nsubs = size(mats,3);
nmods = numel(modules_true);

for ii = 1:nsubs
    mat = mats(:,:,ii);
    mat = weight_conversion(mat,'autofix');
    [~,Q_true(ii),Qc_true(:,ii)] = CR7T_community_louvain_Q_initial(mat,gamma,C,modules_true,B);
end

%% Permutations - 
% Shuffle the random number generator to get a unique set of permutations
rng('shuffle');
% Preallocate results matrix
Qc_perm = zeros(nmods,nperms,nsubs);
% Loop through permutations
% Loop through subjects
for ii = 1:nsubs
    disp(num2str(ii));
    mat = mats(:,:,ii);
    mat = weight_conversion(mat,'autofix');
    parfor jj = 1:nperms
        %disp(num2str(ii));
        perm_vec = randperm(length(C));
        C_perm = C(perm_vec);
        [~,Q_temp(jj),Qc_temp(:,jj)] = CR7T_community_louvain_Q_initial(mat,gamma,C_perm,modules_true,B);
    end 
    Qc_perm(:,:,ii) = Qc_temp;
    Q_perm(:,ii) = Q_temp;
end
    
%% Get significance matrix
for ii = 1:nmods
    for jj = 1:nsubs
        trueval = Qc_true(ii,jj);
        dist = squeeze(Qc_perm(ii,:,jj));
        pval = invprctile(dist,trueval)/100;
        if pval > 0.5
            pval = 1-pval;
        end
        zval = (trueval - mean(dist))/std(dist);
        % Populate full results matrix
        Qc_sig_pval(ii,jj) = pval;
        Qc_sig_zval(ii,jj) = zval;
    end 
end
Qc_sig_hval = zeros(size(Qc_sig_pval));
Qc_sig_hval(Qc_sig_pval < 0.01) = 1;

%% Plot results
% Boxplot
figure('Position',[100,100,1000,500]);
boxplot(Qc_sig_zval'); grid on;
ax = gca;
% Title
ax.Title.String = ['Modularity Contribution (Qc) - Nonparametric Z-score - in ' cond_data ' data - modules from ' cond_modules];
ax.Title.Interpreter = 'none';
% X Ticks
ax.XTick = 1:numel(modules_true);
gcells = num2cell(modules_true);
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 12;
% Axis labels
xL = xlabel('Module #');
xL.FontSize = 12;
ax.YLabel.String = 'Z-score';

% Mean Z
figure('Position',[100,100,1000,500]);
plot(mean(Qc_sig_zval,2),'--.','MarkerSize',30); grid on;
title(['Modularity Contribution (Qc) - Mean nonparametric Z-score across subjects - in ' cond_data ' data - modules from ' cond_modules],'Interpreter','none');
ax = gca;
ax.XTick = 1:numel(modules_true);
gcells = num2cell(modules_true);
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 12;
xL = xlabel('Module #');
xL.FontSize = 12;

% Percentage of subjects significant
figure('Position',[100,100,1000,500]);
plot(mean(Qc_sig_hval,2),'--.','MarkerSize',30); grid on;
title(['Modularity Contribution (Qc) - Percentage of subjects with P < 0.01- in ' cond_data ' data - modules from ' cond_modules],'Interpreter','none');
ax = gca;
ax.XTick = 1:numel(modules_true);
gcells = num2cell(modules_true);
ax.XTickLabel = gcells(1:end);
ax.XTickLabelRotation = 90;
ax.XAxis.FontSize = 12;
xL = xlabel('Module #');
xL.FontSize = 12;


%% Save workspace
cd('Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\module_significance')
save(['workspace_Qc_MODULESfrom_' cond_data '_DATAfrom_' cond_modules '.mat']);

%% Distribution plotting
% n = 5; m = 2;
% plot_pairperm_distribution(Qc_true(m,n),squeeze(Qc_perm(m,:,n)),'Test',2,'Mod2');


