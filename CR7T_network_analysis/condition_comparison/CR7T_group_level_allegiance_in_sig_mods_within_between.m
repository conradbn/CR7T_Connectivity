%% CR7T - Comparison of within and between module allegiance in group level matrices 
% This script performs analysis and plotting regarding values in the
% group-level allegiance matrices. Specifically, this focuses on comparison
% of the allegiance values within and between each module, with modules
% defined seperately for the Symbolic and Nonsymbolic data. The "true" data
% is compared with the "alt" data (i.e. data from the other condition),
% including the difference in overall mean allegiance as well as paired
% t-tests on values from all node pairs within each module. Non-parametric
% z-scores are also computed for these values using pair-permuted data as
% the basis for null distributions.

purge
% Load allegiance matrices
load('workspace_consensus_clustering_gamma_2pt45','Pa_group');
% Load allegiance matrices from permutations
load('workspace_group_consensus_pair_permutations_50000iters_03_27_05_13_22.mat','Pa_group_perm','num_perms');
%load('workspace_group_consensus_pair_permutations_100iters_03_26_22_34_37.mat','Pa_group_perm','num_perms');

% Set Module Colors/Names/IDs
module_id_rgb = {'DMN',2,1,0,0;... % DMN, red
                 'FPN',6,1,1,0;... % FPN, yellow
                 'SMN',10,0,1,1;... % Sensorimotor, cyan
                 'COSN',12,0.8,0,1;... % Cingulo/Opercular/Salience, purple
                 'Aud',16,1,0.7,1;... % STC, pink
                 'CdVNAc',20,1,0.5,0;... % ventral Caudate, orange
                 'Vis',21,0,0,1;... % Visual, blue
                 'Hipp',22,0,0.5,0.5;... % Hippocampus, teal
                 'DAN',24,0,1,0;... % Dorsal Attention (ITG/MTG),green
                 'Prec',27,0.5,0.5,1;... % Precuneus, pale blue
                 'BG',34,0.3,0.1,0;... % Basal ganglia, dark brown
                 'CdD',35,.8,.33,0;... % Caudate, burnt orange
                 'Thal',36,0.5,0.25,0.1;...% Thalumus, light brown
                 'Dropout',256,0.1,0.1,0.1}; % Signal dropout regions

%% ----- SYMBOLIC ---------------------------------------------------------
% Clear everything except the original variables
clearvars -except Pa_group_perm num_perms Pa_group module_id_rgb
close all
% Load variables
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
% Set the distribution plot title string
fig_string = 'Sym > Nonsym Allegiance - Sym Modules'; 

% Run the function
[sym_t,sym_z,sym_pval,sym_z_mean,sym_pval_mean,sym_bf10,sym_perm_tstats] = run_module_allegiance_statistics(Pa_true,Pa_alt,Pa_group_perm,num_perms,sig_modules,C_new,mnames,fig_string);

sym_modules_true = modules_true;
sym_sig_modules = sig_modules;
sym_mnames = mnames; 
sym_mcolors = mcolors;
sym_Qc_sig_pval = Qc_sig_pval;
sym_Qc_sig_zval = Qc_sig_zval;
sym_C_new = C_new;
sym_fig_string = fig_string;

%% ----- NONSYMBOLIC ------------------------------------------------------
% Clear everything except the original variables
clearvars -except Pa_group_perm num_perms Pa_group module_id_rgb sym_* 
close all
% Load variables
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
% Set the distribution plot title string
fig_string = 'Nonsym > Sym Allegiance - Nonsym Modules';
[nonsym_t,nonsym_z,nonsym_pval,nonsym_z_mean,nonsym_pval_mean,nonsym_bf10,nonsym_perm_tstats]= run_module_allegiance_statistics(Pa_true,Pa_alt,Pa_group_perm,num_perms,sig_modules,C_new,mnames,fig_string);

nonsym_modules_true = modules_true;
nonsym_sig_modules = sig_modules;
nonsym_mnames = mnames; 
nonsym_mcolors = mcolors;
nonsym_Qc_sig_pval = Qc_sig_pval;
nonsym_Qc_sig_zval = Qc_sig_zval;
nonsym_C_new = C_new;
nonsym_fig_string = fig_string;

%% Plot the final results 
plot_module_allegiance_statistics(sym_t,sym_z,sym_pval,sym_z_mean,sym_pval_mean,...
                                  sym_bf10,sym_mnames,sym_Qc_sig_zval,sym_modules_true,...
                                  sym_sig_modules,sym_mcolors,num_perms,sym_fig_string);
plot_module_allegiance_statistics(nonsym_t,nonsym_z,nonsym_pval,nonsym_z_mean,nonsym_pval_mean,...
                                  nonsym_bf10,nonsym_mnames,nonsym_Qc_sig_zval,nonsym_modules_true,...
                                  nonsym_sig_modules,nonsym_mcolors,num_perms,nonsym_fig_string);

%% Function
function [t,z,pval,z_mean,pval_mean,bf10,perm_tstats] = run_module_allegiance_statistics(Pa_true,Pa_alt,Pa_group_perm,num_perms,sig_modules,C_new,mnames,fig_string);
    % Calculate the average allegiance within each module for each
    % condition and also get the reduced allegiance matrix values
    for ii = 1:numel(sig_modules)
        for pp = 1:numel(sig_modules)
            inds1 = C_new == sig_modules(ii);
            inds2 = C_new == sig_modules(pp);
            if ii == pp % Within module
                % True matrix (get
                p = Pa_true(inds1,inds2);
                upper = logical(triu(ones(size(p)),1));
                Pa_true_within_btw(ii,pp) = mean(p(upper));
                Pa_true_within_btw_all{ii,pp} = p(upper);
                % Alternative matrix    
                p = Pa_alt(inds1,inds2);
                upper = logical(triu(ones(size(p)),1));
                Pa_alt_within_btw(ii,pp) = mean(p(upper));
                Pa_alt_within_btw_all{ii,pp} = p(upper);
            elseif ii ~= pp % Between module
                % True matrix
                p = Pa_true(inds1,inds2);
                Pa_true_within_btw(ii,pp) = mean(p(:));
                Pa_true_within_btw_all{ii,pp} = p(:);
                % Alternative matrix    
                p = Pa_alt(inds1,inds2);
                upper = logical(triu(ones(size(p)),1));
                Pa_alt_within_btw(ii,pp) = mean(p(:));
                Pa_alt_within_btw_all{ii,pp} = p(:);
            end
        end
    end

    % Construct null distributions from permuted data 
    for ii = 1:numel(sig_modules)
        for pp = 1:numel(sig_modules)
            disp(['Module #' num2str(ii) ' x Module #' num2str(pp)]);
            inds1 = C_new == sig_modules(ii);
            inds2 = C_new == sig_modules(pp);
            p1 = Pa_group_perm(inds1,inds2,:,1);
            p2 = Pa_group_perm(inds1,inds2,:,2);
            parfor jj = 1:num_perms
                if ii == pp % Within module
                    % Matrix 1
                    p = p1(:,:,jj);
                    upper = logical(triu(ones(size(p)),1));
                    Pa_perm1(jj) = mean(p(upper));
                    Pa_perm1_all = p(upper);
                    % Matrix 2    
                    p = p2(:,:,jj);
                    upper = logical(triu(ones(size(p)),1));
                    Pa_perm2(jj) = mean(p(upper));
                    Pa_perm2_all = p(upper);
                elseif ii ~= pp % Between module
                    % Matrix 1
                    p = p1(:,:,jj);
                    Pa_perm1(jj) = mean(p(:));
                    Pa_perm1_all = p(:);
                    % Matrix 2
                    p = p2(:,:,jj);
                    Pa_perm2(jj) = mean(p(:));
                    Pa_perm2_all = p(:);
                end
                % Run paired t-tests
                [~,p,~,stat] = ttest(Pa_perm1_all,Pa_perm2_all);
                tstat(jj) =  single(floor(stat.tstat*10000)/10000); % ROUND TO ACCOUNT FOR PRECISION ERRORS
                
%                 if isequal(Pa_true_within_btw_all{ii,pp},Pa_alt_within_btw_all{ii,pp})
%                     tstat(jj) = 0;
%                 elseif length(Pa_true_within_btw_all{ii}) == 1
%                     tstat(jj) = 0;
%                 else
%                     [~,p,~,stat]=ttest(Pa_perm1_all,Pa_perm2_all);
%                     tstat(jj) = stat.tstat;
%                 end
            end
            perm_tstats(ii,pp,:) = tstat; 
            Pa_perm1_within_btw(ii,pp,:) = Pa_perm1;
            Pa_perm2_within_btw(ii,pp,:) = Pa_perm2;
        end
    end  

    % Get difference matrix
    Pa_perm_within_btw_diff = Pa_perm1_within_btw - Pa_perm2_within_btw;

    % Calculate z-score of t-stat for each modules
    for ii = 1:numel(sig_modules)
        for pp = 1:numel(sig_modules)
            if isequal(Pa_true_within_btw_all{ii,pp},Pa_alt_within_btw_all{ii,pp}) % Don't run ttest on the same data as it gives NaN
                t(ii,pp) = 0;
            elseif length(Pa_true_within_btw_all{ii,pp}) == 1
                t(ii,pp) = 0;
            else
                [bf10(ii,pp),pbf] = bf.ttest(Pa_true_within_btw_all{ii,pp},Pa_alt_within_btw_all{ii,pp});
                [~,p(ii,pp),~,stat] = ttest(Pa_true_within_btw_all{ii,pp},Pa_alt_within_btw_all{ii,pp});
                t(ii,pp) = single(floor(stat.tstat*10000)/10000); %ROUND TO ACCOUNT FOR PRECISION ERRORS
                z(ii,pp) = (stat.tstat-mean(perm_tstats(ii,pp,:)))/std(perm_tstats(ii,pp,:));                   
                pval(ii,pp) = sum(squeeze(perm_tstats(ii,pp,:) <= t(ii,pp)))/num_perms;% invprctile(perm_tstats(ii,pp,:),stat.tstat); INVPRCTILE DOESNT WORK WITH REPEATED VALUES
                if pval(ii,pp) == 1
                    pval(ii,pp) = sum(squeeze(perm_tstats(ii,pp,:) < t(ii,pp)))/num_perms; % Account for repeated max value
                end
                % Plot true value on distribution
                fig = figure; 
                plot_pairperm_distribution(stat.tstat,squeeze(perm_tstats(ii,pp,:)),[fig_string ' - ' mnames{ii} ' x ' mnames{pp}],4,0); xlabel('T-stat');
                ax = gca;
                ax.YLim = [0 .15];
                out = strrep(ax.Title.String,' ','_');
                out = strrep(out,'>','gr');
                print(fig,out,'-dpng');
            end
        % Get the difference in mean allegiance (for case where only two
        % regions/one pairwise value)
        z_mean(ii,pp) = (Pa_true_within_btw(ii,pp) - Pa_alt_within_btw(ii,pp)- mean(Pa_perm_within_btw_diff(ii,pp,:)))...
                        /std(Pa_perm_within_btw_diff(ii,pp,:));
        pval_mean(ii,pp) = sum(squeeze(Pa_perm_within_btw_diff(ii,pp,:) < ...
                           (Pa_true_within_btw(ii,pp)-Pa_alt_within_btw(ii,pp))))/num_perms;%invprctile(Pa_perm_within_btw_diff(ii,pp,:),Pa_true_within_btw(ii,pp)-Pa_alt_within_btw(ii,pp));
        end
    end
end

%% Plotting
function [] = plot_module_allegiance_statistics(t,z,pval,z_mean,pval_mean,bf10,mnames,Qc_sig_zval,modules_true,sig_modules,mcolors,num_perms,fig_string)

    % Get sorted order of data for plotting (sorted based on subject-level mean Q*c z-scores, as in initial figure)
    [~,msort] = sort(mean(Qc_sig_zval,2),'descend');
    msort = modules_true(msort);
    [~,msort] = ismember(msort,sig_modules); % Find the location of significant modules in sorted module list
    msort(msort == 0) = []; % Remove 0s, i.e. nonsignificant modules

    % Actually sort the data now
    t = t(msort,msort);
    z = z(msort,msort);
    bf10 = bf10(msort,msort);
    pval = pval(msort,msort);
    z_mean = z_mean(msort,msort);
    pval_mean = pval_mean(msort,msort);
    mnames = mnames(msort);
    mcolors = mcolors(msort,:);

    for ii = 1:numel(mnames)
        if strcmp(mnames(ii),'CaudD')
            z(ii,ii) = z_mean(ii,ii);
            pval(ii,ii) = pval_mean(ii,ii);
        end
    end
    close all

    RdBu = cbrewer('div', 'RdBu', 256);
    RdBu = flip(RdBu);
    % figure;
    % imagesc(z);
    % ax = gca;
    % ax.XTickLabel = mnames;
    % ax.YTickLabel = mnames;
    % colormap(RdBu);
    
    %DONT NEED TO DIVIDE BY 100 SINCE NOT USiNG INVPRCTILE
        for ii = 1:size(pval,1)
            for jj = 1:size(pval,2)
                if pval(ii,jj) < 0.5
                    pval2(ii,jj) = pval(ii,jj);
                else
                    pval2(ii,jj) = 1-pval(ii,jj);
                end
                if pval2(ii,jj) == 0
                    pval2(ii,jj) = 1;
                end
            end
        end

    % Plot the diagonal and below (z-scores)
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
    ax.FontSize = 20;
    hold on
    %axis square
    ax.Title.String = fig_string;
    
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
    % Bonferoni
    %pval2_lower_corrected = pval2_lower * sum(~isnan(z_lower(:)));
%     
%     % P<0.05 (two-tailed)
%     [row,col] = find(pval2_lower<0.025);
%     for ii = 1:numel(row)
%         if t(row(ii),col(ii)) == 0
%             continue
%         else
%             plot(col(ii),row(ii),'k.','markers', 30);
%         end
%     end
%     % P<0.01 (two-tailed)
%     [row,col] = find(pval2_lower<0.005);
%     for ii = 1:numel(row)
%         if t(row(ii),col(ii)) == 0
%             continue
%         else
%             plot(col(ii),row(ii),'k.','markers', 45);
%         end
%     end
%     % FDR P<0.05
%     [row,col] = find(pval2_lower_corrected <0.05);
%     for ii = 1:numel(row)
%         if t(row(ii),col(ii)) == 0
%             continue
%         else
%             plot(col(ii),row(ii),'k.','markers', 65);
%         end
%     end

    % Title and print png
    out = strrep(ax.Title.String,' ','_');
    out = strrep(out,'>','gr');
    print(fig,out,'-dpng');

    pval2_diag = diag(pval2);
    [~, ~, ~, adj_p] = fdr_bh(pval2_diag,.05,'dep','yes');
end