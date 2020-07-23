%% Zscore connectivity matrix permutations
% 1/21/20 
% Get pair-permuted subject-level connectivity matrices. Can use
% output for testing differences between Symbolic and Nonsymbolic
% conditions in true wholebrain/module-level/region-level metrics

tic
purge;
% Shuffle the random number generator to get a unique set of permutations
rng('shuffle');

% Start up the parallel pool of CPU "workers"
%parpool

% Set number of permutations to run
num_perms = 50000; 

% Load subject connectivity matrices
load('workspace_zscore_mats_33subs.mat','zscore_mats');

% Set connectivity matrices for permutations
zscore_perms(:,:,:,1) = single(zscore_mats.BN202_Compare_Symbolic);
zscore_perms(:,:,:,2) = single(zscore_mats.BN202_Compare_Nonsymbolic);
num_subs = size(zscore_perms,3);

%% Get community information
% Set Module Colors/Names/IDs
module_id_rgb = {'DMN',2,1,0,0;... % DMN, red
                 'FPN',6,1,1,0;... % FPN, yellow
                 'SMN',10,0,1,1;... % Sensorimotor, cyan
                 'SN',12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
                 'STC',16,1,0.7,1;... % STC, pink
                 'CaudV',20,1,0.5,0;... % ventral Caudate, orange
                 'Vis',21,0,0,1;... % Visual, blue
                 'Hipp',22,0,0.5,0.5;... % Hippocampus, teal
                 'DAN',24,0,1,0;... % Dorsal Attention (ITG/MTG),green
                 'Prec',27,0.5,0.5,1;... % Precuneus, pale blue
                 'BG',34,0.3,0.1,0;... % Basal ganglia, dark brown
                 'CaudD',35,.8,.33,0;... % Caudate, burnt orange
                 'Thal',36,0.5,0.25,0.1;...% Thalumus, light brown
                 'Dropout',256,0.1,0.1,0.1}; % Signal dropout regions
% Load variables
load('workspace_Qc_significance_symbolic.mat','Qc_sig_pval','Qc_sig_zval','modules_true','C');
% Get significant modules based on p-value (since some of the permuted
% distributions are not normal,i.e. have long positive tail)
sig_modules = modules_true(median(Qc_sig_pval,2) < 0.01);
C_new = C .* ismember(C,sig_modules);

% Set the module colors and names from the defined cell array
ix = 1;
for ii = 1:length(modules_true)
    if ismember(modules_true(ii),cell2mat(module_id_rgb(:,2))) && ismember(modules_true(ii),sig_modules)
        ind = find(cell2mat(module_id_rgb(:,2)) == modules_true(ii));
        mcolors(ix,:) = cell2mat(module_id_rgb(ind,3:5));
        mnames(ix) = module_id_rgb(ind,1);
        ix = ix + 1;
    end
end

%% Compute true pairwise connectivity differences within/between communities
for ii = 1:numel(sig_modules)
    for jj = 1:numel(sig_modules)
        inds1 = C_new == sig_modules(ii);
        inds2 = C_new == sig_modules(jj);
        a = zscore_mats.BN202_Compare_Symbolic(inds1,inds2,:);
        amed = median(a,3);
        b = zscore_mats.BN202_Compare_Nonsymbolic(inds1,inds2,:);
        bmed = median(b,3);
        if ii == jj
            upper = logical(triu(ones(size(amed)),1));
            amed = amed(upper);
            bmed = bmed(upper);
        end
        [~,pval,~,stat] = ttest(amed(:),bmed(:));
        t_true(ii,jj) = stat.tstat;
        pval_true(ii,jj) = pval; 
    end
end

%% Construct null distribution using permutation
disp(['Working on ' num2str(num_perms) ' permutations']);
% Progress bar
%parfor_progress(num_perms);

waypoints = 0:100:num_perms;

parfor p = 1:num_perms
    if ismember(p,waypoints)
        disp(['On permutation # ' num2str(p) ' of ' num2str(num_perms) ' ...']);
    end
    % Get randomly relabeled set of matrices
    inds_cond1 = randi([1 2],num_subs,1);
    inds_cond2 = inds_cond1-1;
    inds_cond2(inds_cond2 == 0) = 2;
    z1 = [];
    z2 = [];
    for s = 1:num_subs
        z1(:,:,s) = zscore_perms(:,:,s,inds_cond1(s));
        z2(:,:,s) = zscore_perms(:,:,s,inds_cond2(s));
    end
    % Loop through communities
    t_tmp = [];
    for ii = 1:numel(sig_modules)
        for jj = 1:numel(sig_modules)
            inds1 = C_new == sig_modules(ii);
            inds2 = C_new == sig_modules(jj);
            a = z1(inds1,inds2,:);
            amed = median(a,3);
            b = z2(inds1,inds2,:);
            bmed = median(b,3);
            if ii == jj
                upper = logical(triu(ones(size(amed)),1));
                amed = amed(upper);
                bmed = bmed(upper);
            end
            % Compute ttest on median connectivity vectors
            [~,~,~,stat] = ttest(amed(:),bmed(:));
            t_tmp(ii,jj) = stat.tstat;
        end
    end
    t_perm(:,:,p) = t_tmp;
    % Print current progress
    %parfor_progress;
end

%% Calculate the z-score and non-parametric p-value from null distribution
for ii = 1:size(t_true,1)
    for jj = 1:size(t_true,2)
        t = t_true(ii,jj); 
        null_dist = squeeze(t_perm(ii,jj,:));
        zscore_true(ii,jj) = (t-mean(null_dist))/std(null_dist);
        npval_true(ii,jj) = sum(null_dist < t)/num_perms;
    end
end

%% Sort the data for manuscript display
% Get sorted order of data for plotting (sorted based on subject-level mean Q*c z-scores, as in initial figure)
[~,msort] = sort(mean(Qc_sig_zval,2),'descend');
msort = modules_true(msort);
[~,msort] = ismember(msort,sig_modules); % Find the location of significant modules in sorted module list
msort(msort == 0) = []; % Remove 0s, i.e. nonsignificant modules

zscore_true = zscore_true(msort,msort);
npval_true = npval_true(msort,msort);
mnames = mnames(msort);
mcolors = mcolors(msort,:);

%% Plot results matrices
RdBu = cbrewer('div', 'RdBu', 256);
RdBu = flip(RdBu);

% Fix Pvals
for ii = 1:size(npval_true,1)
    for jj = 1:size(npval_true,2)
        if npval_true(ii,jj) < 0.5
            pval2(ii,jj) = npval_true(ii,jj);
        else
            pval2(ii,jj) = 1-npval_true(ii,jj);
        end
%         if pval2(ii,jj) == 0
%             pval2(ii,jj) = 1;
%         end
    end
end

% Plot the diagonal and below (z-scores)
z = zscore_true;
upper = logical(triu(ones(size(z)),1));
z_lower = z;
z_lower(upper) = NaN;
fig = figure('Position',[100,100,800,600]);
im = imagesc(z_lower,'AlphaData',double(~isnan(z_lower)));
ax = gca;
ax.Color = [0.8,0.8,0.8];
ax.XTick = 1:numel(mnames);
ax.YTick = 1:numel(mnames);
ax.XTickLabel = mnames;
ax.YTickLabel = mnames;
colormap(RdBu); colorbar
grid on
ax.CLim = [-2.5 2.5];
ax.FontSize = 15;
hold on
ax.Title.String = 'Sym > Nonsym Connectivity - Group Median - in Sym Modules';

% Plot boxlines
% Get module box lines
X = [];
Y = [];
for i =1:size(z,1)
    mn = min(i) - 0.5;
    mx = max(i) + 0.5;
    x = [mn mn mx mx mn NaN];
    y = [mn mx mx mn mn NaN];
    X = [X, x]; 
    Y = [Y, y];
end
plot(X,Y,'k','Linewidth',2);

% Plot significance indicators
pval2_lower = pval2;
pval2_lower(upper) = 1;

% Get FDR adjbusted pvals
pval2_lower_vec = pval2_lower(:);
inds_use = find(pval2_lower_vec ~=1);
[~, ~, ~, adj_p] = fdr_bh(pval2_lower_vec(inds_use),.05,'pdep','yes');
pval2_lower_vec(inds_use) = adj_p;
pval2_lower_corrected = reshape(pval2_lower_vec,[size(pval2,1) size(pval2,2)]);

% P<0.05 (two-tailed)

[row,col] = find(pval2_lower<0.025);
for ii = 1:numel(row)
    if t_true(row(ii),col(ii)) == 0
        continue
    else
        plot(col(ii),row(ii),'k.','markers', 25);
    end
end
% P<0.01 (two-tailed)
[row,col] = find(pval2_lower<0.005);
for ii = 1:numel(row)
    if t_true(row(ii),col(ii)) == 0
        continue
    else
        plot(col(ii),row(ii),'kd','markers', 12,'MarkerFaceColor','k');
    end
end
% FDR P<0.05
[row,col] = find(pval2_lower_corrected <0.05);
for ii = 1:numel(row)
    if t_true(row(ii),col(ii)) == 0
        continue
    else
        plot(col(ii),row(ii),'kp','markers', 15,'MarkerFaceColor','k');
    end
end

% Title and print png
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
%print(fig,out,'-dpng');

pval2_diag = diag(pval2);
[~, ~, ~, adj_p] = fdr_bh(pval2_diag,.05,'dep','yes');

%% Save data
save(['Z:\CR7T_Connectivity\matfiles\CR7T_WORKSPACES\'...
      'workspace_median_connectivity_comparison_sym_comms_' num2str(num_perms) 'iters_' datestr(now,'mm_DD_HH_MM_SS')],'zscore_true','npval_true','mnames','-v7.3');
% Report total compute time
toc








