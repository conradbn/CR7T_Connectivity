%% CR7T - Comparison of within-module allegiance in group level matrices 
% This script performs analysis and plotting regarding values in the
% group-level allegiance matrices. Specifically, this focuses on comparison of the
% allegiance values within each module, with modules defined seperately for
% the Symbolic and Nonsymbolic data. The "true" data is compared with the
% "alt" data (i.e. data from the other condition), including the difference
% in overall mean allegiance as well as paired t-tests on values from all node
% pairs within each module. Non-parametric z-scores are also computed for these 
% values using pair-permuted data as the basis for null distributions.

purge
% Load allegiance matrices
load('workspace_consensus_clustering_gamma_2pt45','Pa_group');
% Load allegiance matrices from permutations
load('workspace_group_consensus_pair_permutations_50000iters_03_27_05_13_22.mat','Pa_group_perm','num_perms');

% Set Module Colors/Names/IDs
module_id_rgb = {'DMN',2,1,0,0;... % DMN, red
                 'FPN',6,1,1,0;... % FPN, yellow
                 'SMN',10,0,1,1;... % Sensorimotor, cyan
                 'Sal',12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
                 'STS',16,1,0.7,1;... % STS, dark gray
                 'CaudV',20,1,0.5,0;... % ventral Caudate, orange
                 'Vis',21,0,0,1;... % Visual, blue
                 'Hipp',22,0,0.5,0.5;... % Hippocampus, teal
                 'DAN',24,0,1,0;... % Dorsal Attention (ITG/MTG),green
                 'Prec',27,0.5,0.5,1;... % Precuneus, pale blue
                 'BG',34,0.3,0.1,0;... % Basal ganglia, dark brown
                 'CaudD',35,.8,.33,0;... % Caudate, burnt orange
                 'Thal',36,0.5,0.25,0.1;...% Thalumus, light brown
                 'Dropout',256,0.1,0.1,0.1}; % Signal dropout regions

%% Symbolic ---------------------------------------------------------------------
load('workspace_Qc_significance_symbolic.mat','Qc_sig_pval','Qc_sig_zval','modules_true','C');

% Get significant modules based on p-value (since some of the permuted distributions are not normal,i.e. have long positive tail)
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

% Set the allegiance matrices to work with;
Pa_true = Pa_group.BN202_Compare_Symbolic; 
Pa_alt = Pa_group.BN202_Compare_Nonsymbolic;

% Calculate the average connectivity within each module for each
% condition
for ii = 1:numel(sig_modules)
    inds = C_new == sig_modules(ii);
    % True matrix
    p = Pa_true(inds,inds);
    upper = logical(triu(ones(size(p)),1));
    Pa_true_within(ii) = mean(p(upper));
    Pa_true_within_all{ii} = p(upper);
    % Alternative matrix    
    p = Pa_alt(inds,inds);
    upper = logical(triu(ones(size(p)),1));
    Pa_alt_within(ii) = mean(p(upper));
    Pa_alt_within_all{ii} = p(upper);
end

% Construct null distributions from permuted data 
for ii = 1:numel(sig_modules)
    disp(['Module #' num2str(ii)]);
    inds = C_new == sig_modules(ii);
    for jj = 1:num_perms
        % Matrix 1
        p = Pa_group_perm(inds,inds,jj,1);
        upper = logical(triu(ones(size(p)),1));
        Pa_perm1_within(ii,jj) = mean(p(upper));
        Pa_perm1_within_all = p(upper);
        % Matrix 2    
        p = Pa_group_perm(inds,inds,jj,2);
        upper = logical(triu(ones(size(p)),1));
        Pa_perm2_within(ii,jj) = mean(p(upper));
        Pa_perm2_within_all = p(upper);
        % Run paired t-tests
        if isequal(Pa_true_within_all{ii},Pa_alt_within_all{ii})%isequal(Pa_perm1_within_all{ii},Pa_perm2_within_all{ii})
            t_perm(ii,jj) = 0;
        elseif length(Pa_true_within_all{ii}) == 1%length(Pa_perm1_within_all{ii}) == 1
            t_perm(ii,jj) = 0;
        else
            [~,p,~,stat]=ttest(Pa_perm1_within_all,Pa_perm2_within_all);
            t_perm(ii,jj) = stat.tstat;
        end
    end
end  

Pa_perm_within_diff = Pa_perm1_within - Pa_perm2_within;

% Calculate z-score of t-stat for each modules
for ii = 1:numel(sig_modules)
    if isequal(Pa_true_within_all{ii},Pa_alt_within_all{ii})
        t(ii) = 0;
    elseif length(Pa_true_within_all{ii}) == 1
        t(ii) = 0;
    else
        [bf10(ii),pbf] = bf.ttest(Pa_true_within_all{ii},Pa_alt_within_all{ii});
        [~,p(ii),~,stat]=ttest(Pa_true_within_all{ii},Pa_alt_within_all{ii});
%         [psr(ii),~,statsr]=signrank(Pa_true_within_all{ii},Pa_alt_within_all{ii});
%         zvalsr(ii) = statsr.zval;
        t(ii) = stat.tstat;
        z(ii) = (stat.tstat-mean(t_perm(ii,:)))/std(t_perm(ii,:));
        pval(ii) = invprctile(t_perm(ii,:),stat.tstat);
        % Plot true value on distribution
        fig = figure; 
        plot_pairperm_distribution(stat.tstat,t_perm(ii,:),['Sym > Nonsym Allegiance - Symbolic Module - ' mnames{ii}],4,0); xlabel('T-stat');
        ax = gca;
        ax.YLim = [0 .15];
        out = strrep(ax.Title.String,' ','_');
        out = strrep(out,'>','gr');
        print(fig,out,'-dpng');
    end
    z_mean(ii) = (Pa_true_within(ii)-Pa_alt_within(ii)-mean(Pa_perm_within_diff(ii,:)))/std(Pa_perm_within_diff(ii,:));
    pval_mean(ii) = invprctile(Pa_perm_within_diff(ii,:),Pa_true_within(ii)-Pa_alt_within(ii));
end

% Get sorted order of data for plotting (sorted based on subject-level mean Q*c z-scores, as in initial figure)
[~,msort] = sort(mean(Qc_sig_zval,2),'descend');
msort = modules_true(msort);
[~,msort] = ismember(msort,sig_modules); % Find the location of significant modules in sorted module list
msort(msort == 0) = []; % Remove 0s, i.e. nonsignificant modules

% Actually sort the data now
t = t(msort);
z = z(msort);
bf10 = bf10(msort);
z_mean = z_mean(msort);
mnames = mnames(msort);
mcolors = mcolors(msort,:);

% Plot T-stat of paired T-test
fig = figure;
br = bar(t); grid on; 
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
ax.Title.String = 'T-stat - Sym > Nonsym - Paired t-test - Allegiance';
ax.YLim = [-6 14];
ax.YLabel.String = 'T-stat';
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 300;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Plot nonparametric Z-score of t-stat of paired T-test
fig = figure('Position',[100,100,1400,800]);
br = bar(z); grid on;
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
%ax.Title.String = 'Z-score of T-stat - Sym > Nonsym - Paired t-test - Allegiance';
ax.Title.String = 'Region-to-Region Integration - Symbolic > Nonsymbolic';
ax.YLim = [-1.5 3];
ax.YLabel.String = 'Z-score of T-stat (paired test)';
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 280;
ax.FontSize = 25;
ax.GridAlpha = .5;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Plot nonparametric Z-score of mean difference 
fig = figure;
br = bar(z_mean); grid on;
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
ax.Title.String = 'Z-score of Difference in Mean Allegiance - Sym > Nonsym';
ax.YLim = [-1.5 3];
ax.YLabel.String = 'Z-score of Difference';
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 300;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Plot Bayes Factor of paired T-test
fig = figure;
br = bar(bf10); grid on;
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
ax.Title.String = 'Bayes Factor (BF10) - Sym > Nonsym - Paired t-test - Allegiance';
ax.YScale = 'log';
ax.YLabel.String = 'BF10';
%ax.YLim = [-1.5 3];
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 300;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Barcharts of the true and alternative allegiance values from each module
% First create full matrix for boxplotting, padding smaller modules with NaNs
Pa_true_within_all_sort = Pa_true_within_all(msort);
Pa_alt_within_all_sort = Pa_alt_within_all(msort);
Pa_true_alt = [];
for ii = 1:numel(sig_modules)
    Ptrue = Pa_true_within_all_sort{ii};
    Palt = Pa_alt_within_all_sort{ii};
    num_nodes = length(Ptrue);
    if ii~=1 && num_nodes~=length(Pa_true_alt)
        Ptrue = [Ptrue; nan(length(Pa_true_alt)-num_nodes,1)];
        Palt = [Palt; nan(length(Pa_true_alt)-num_nodes,1)];
    end
    Pa_true_alt = [Pa_true_alt,Ptrue,Palt];
end

% Custom bar chart
data = Pa_true_alt;
data_mean = nanmean(data);
for ii = 1:size(data,2)
    std_err_mean(1,ii) = nanstd(data(:,ii))/sqrt(numel(~isnan(data(:,ii))));
end
figure('Position',[100,100,1400,600]);
br = bar(data_mean);
br.BarWidth = 0.9;
br.FaceColor = 'flat';
mcolors2 = repelem(mcolors,2,1);
for ii = 1:length(mcolors2)
    br.CData(ii,:)= mcolors2(ii,:);
end
% Axis properties
ax = gca;
ax.YLim = [0 1];
ax.YLabel.String = 'Mean Allegiance';
ax.XTick = 1:size(data,2);
ax.XTickLabel = repelem(mnames,2,1);
ax.XTickLabelRotation = 280;
ax.Title.String = 'Mean Allegiance in Symbolic Modules - Symbolic vs Nonsymbolic';
ax.FontSize = 25;
ax.GridAlpha = .5;
% Error bars
hold on
er = errorbar(1:size(data,2),data_mean,std_err_mean,std_err_mean);    
er.Color = 'k';
er.LineWidth = 3;
er.LineStyle = 'none';  
grid on
hold off
print(gcf,strrep(ax.Title.String,' ','_'),'-dpng');

%% Clear everything except the original variables
clearvars -except Pa_group_perm num_perms Pa_group module_id_rgb
close all

%% Nonsymbolic ---------------------------------------------------------------------
load('workspace_Qc_significance_nonsymbolic.mat','Qc_sig_pval','Qc_sig_zval','modules_true','C');

% Get significant modules based on p-value (since some of the permuted distributions are not normal,i.e. have long positive tail)
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

% Set the allegiance matrices to work with;
Pa_true = Pa_group.BN202_Compare_Nonsymbolic; 
Pa_alt = Pa_group.BN202_Compare_Symbolic;

% Calculate the average connectivity within each module for each
% condition
for ii = 1:numel(sig_modules)
    inds = C_new == sig_modules(ii);
    % True matrix
    p = Pa_true(inds,inds);
    upper = logical(triu(ones(size(p)),1));
    Pa_true_within(ii) = mean(p(upper));
    Pa_true_within_all{ii} = p(upper);
    % Alternative matrix    
    p = Pa_alt(inds,inds);
    upper = logical(triu(ones(size(p)),1));
    Pa_alt_within(ii) = mean(p(upper));
    Pa_alt_within_all{ii} = p(upper);
end

% Construct null distributions from permuted data 
for ii = 1:numel(sig_modules)
    disp(['Module #' num2str(ii)]);
    inds = C_new == sig_modules(ii);
    for jj = 1:num_perms
        % Matrix 1
        p = Pa_group_perm(inds,inds,jj,1);
        upper = logical(triu(ones(size(p)),1));
        Pa_perm1_within(ii,jj) = mean(p(upper));
        Pa_perm1_within_all = p(upper);
        % Matrix 2    
        p = Pa_group_perm(inds,inds,jj,2);
        upper = logical(triu(ones(size(p)),1));
        Pa_perm2_within(ii,jj) = mean(p(upper));
        Pa_perm2_within_all = p(upper);
        % Run paired t-tests
        if isequal(Pa_true_within_all{ii},Pa_alt_within_all{ii})%isequal(Pa_perm1_within_all{ii},Pa_perm2_within_all{ii})
            t_perm(ii,jj) = 0;
        elseif length(Pa_true_within_all{ii}) == 1%length(Pa_perm1_within_all{ii}) == 1
            t_perm(ii,jj) = 0;
        else
            [~,p,~,stat]=ttest(Pa_perm1_within_all,Pa_perm2_within_all);
            t_perm(ii,jj) = stat.tstat;
        end
    end
end  

Pa_perm_within_diff = Pa_perm1_within - Pa_perm2_within;

% Calculate z-score of t-stat for each modules
for ii = 1:numel(sig_modules)
    if isequal(Pa_true_within_all{ii},Pa_alt_within_all{ii})
        t(ii) = 0;
    elseif length(Pa_true_within_all{ii}) == 1
        t(ii) = 0;
    else
        [bf10(ii),pbf] = bf.ttest(Pa_true_within_all{ii},Pa_alt_within_all{ii});
        [~,p(ii),~,stat]=ttest(Pa_true_within_all{ii},Pa_alt_within_all{ii});
        t(ii) = stat.tstat;
        z(ii) = (stat.tstat-mean(t_perm(ii,:)))/std(t_perm(ii,:));
        pval(ii) = invprctile(t_perm(ii,:),stat.tstat);
        % Plot true value on distribution
        fig = figure; 
        plot_pairperm_distribution(stat.tstat,t_perm(ii,:),['Nonsym > Sym Allegiance - Nonsymbolic Module - ' mnames{ii}],4,0); xlabel('T-stat');
        ax = gca;
        ax.YLim = [0 .15];
        out = strrep(ax.Title.String,' ','_');
        out = strrep(out,'>','gr');
        print(fig,out,'-dpng');
    end
    z_mean(ii) = (Pa_true_within(ii)-Pa_alt_within(ii)-mean(Pa_perm_within_diff(ii,:)))/std(Pa_perm_within_diff(ii,:));
    pval_mean(ii) = invprctile(Pa_perm_within_diff(ii,:),Pa_true_within(ii)-Pa_alt_within(ii));
end

% Get sorted order of data for plotting (sorted based on subject-level mean Q*c z-scores, as in initial figure)
[~,msort] = sort(mean(Qc_sig_zval,2),'descend');
msort = modules_true(msort);
[~,msort] = ismember(msort,sig_modules); % Find the location of significant modules in sorted module list
msort(msort == 0) = []; % Remove 0s, i.e. nonsignificant modules

% Actually sort the data now
t = t(msort);
z = z(msort);
bf10 = bf10(msort);
z_mean = z_mean(msort);
mnames = mnames(msort);
mcolors = mcolors(msort,:);

% Plot T-stat of paired T-test
fig = figure;
br = bar(t); grid on; 
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
ax.Title.String = 'T-stat - Nonsym > Sym - Paired t-test - Allegiance';
ax.YLim = [-6 14];
ax.YLabel.String = 'T-stat';
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 300;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Plot nonparametric Z-score of t-stat of paired T-test
fig = figure('Position',[100,100,1400,800]);
  % Insert CaudD z-score of Difference (since can't do paired test with
  % only one connection, i.e. only two regions)
  z2 = z;
  z2(10) = z_mean(10);

br = bar(z2); grid on;
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
%ax.Title.String = 'Z-score of T-stat - Nonsym > Sym - Paired t-test - Allegiance';
ax.Title.String = 'Region-to-Region Integration - Nonsymbolic > Symbolic';
ax.YLim = [-1.5 3];
ax.YLabel.String = 'Z-score of T-stat (paired test)';
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 280;
ax.FontSize = 25;
ax.GridAlpha = .5;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Plot nonparametric Z-score of mean difference 
fig = figure;
br = bar(z_mean); grid on;
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
ax.Title.String = 'Z-score of Difference in Mean Allegiance - Nonsym > Sym';
ax.YLim = [-1.5 3];
ax.YLabel.String = 'Z-score of Difference';
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 300;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');
% Plot Bayes Factor of paired T-test
fig = figure;
br = bar(bf10); grid on;
br.FaceColor = 'flat';
br.CData = mcolors;
ax = gca;
ax.Title.String = 'Bayes Factor (BF10) - Nonsym > Sym - Paired t-test - Allegiance';
ax.YScale = 'log';
ax.YLabel.String = 'BF10';
%ax.YLim = [-1.5 3];
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 300;
out = strrep(ax.Title.String,' ','_');
out = strrep(out,'>','gr');
print(fig,out,'-dpng');

% Barcharts of the true and alternative allegiance values from each module
% First create full matrix for boxplotting, padding smaller modules with NaNs
Pa_true_within_all_sort = Pa_true_within_all(msort);
Pa_alt_within_all_sort = Pa_alt_within_all(msort);
Pa_true_alt = [];
for ii = 1:numel(sig_modules)
    Ptrue = Pa_true_within_all_sort{ii};
    Palt = Pa_alt_within_all_sort{ii};
    num_nodes = length(Ptrue);
    if ii~=1 && num_nodes~=length(Pa_true_alt)
        Ptrue = [Ptrue; nan(length(Pa_true_alt)-num_nodes,1)];
        Palt = [Palt; nan(length(Pa_true_alt)-num_nodes,1)];
    end
    Pa_true_alt = [Pa_true_alt,Ptrue,Palt];
end

% Custom bar chart
data = Pa_true_alt;
data_mean = nanmean(data);
for ii = 1:size(data,2)
    std_err_mean(1,ii) = nanstd(data(:,ii))/sqrt(numel(~isnan(data(:,ii))));
end
figure('Position',[100,100,1400,600]);
br = bar(data_mean);
br.BarWidth = 0.9;
br.FaceColor = 'flat';
mcolors2 = repelem(mcolors,2,1);
for ii = 1:length(mcolors2)
    br.CData(ii,:)= mcolors2(ii,:);
end
% Axis properties
ax = gca;
ax.YLim = [0 1];
ax.YLabel.String = 'Mean Allegiance';
ax.XTick = 1:size(data,2);
ax.XTickLabel = repelem(mnames,2,1);
ax.XTickLabelRotation = 280;
ax.Title.String = 'Mean Allegiance in Nonsymbolic Modules - Nonsymbolic vs Symbolic';
ax.FontSize = 25;
ax.GridAlpha = .5;
% Error bars
hold on
er = errorbar(1:size(data,2),data_mean,std_err_mean,std_err_mean);    
er.Color = 'k';
er.LineWidth = 3;
er.LineStyle = 'none';  
grid on
hold off
print(gcf,strrep(ax.Title.String,' ','_'),'-dpng');
