%% CR7T Qc significance plotting
purge;
cd('Z:\CR7T_Connectivity\matfiles\CR7T_network_analysis\module_significance');
%% Set module names and colors
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

%% ~~~~~~ SYMBOLIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Load data
load('workspace_Qc_MODULESfrom_BN202_Compare_Symbolic_DATAfrom_BN202_Compare_Symbolic.mat',...
     'Qc_sig_zval','Qc_sig_pval','Qc_sig_hval','Qc_perm','Qc_true','cond_data','cond_modules','modules_true');
 
%% Get significant modules based on p-value (since some of the permuted distributions are not normal,i.e. have long positive tail)
sig_modules = modules_true(median(Qc_sig_pval,2) < 0.05);

%% Sort the modules by mean z-score
[~,mz_sort] = sort(mean(Qc_sig_zval,2),'descend');
modules_true_sort = modules_true(1,mz_sort);
Qc_sig_zval_sort = Qc_sig_zval(mz_sort,:);

%% Plot selected distributions
% for ii = 1:size(Qc_perm,3)
%     d = squeeze(Qc_perm(13,:,ii));
%     t = Qc_true(13,ii);
%     figure;
%     plot_pairperm_distribution(t,d,'',1,'');
% end

%% Plot results
% Custom bar chart
data = Qc_sig_zval_sort';
data_mean = mean(data);
for ii = 1:size(data,2)
    std_err_mean(1,ii) = nanstd(data(:,ii))/sqrt(numel(~isnan(data(:,ii))));
end
figure('Position',[100,100,1400,800]);
br = bar(data_mean);
br.BarWidth = 0.9;
br.FaceColor = 'flat';
% Set the module colors and names from the defined cell array
for ii = 1:length(data_mean)
    if ismember(modules_true_sort(ii),cell2mat(module_id_rgb(:,2))) && ismember(modules_true_sort(ii),sig_modules)
        ind = find(cell2mat(module_id_rgb(:,2)) == modules_true_sort(ii));
        mcolors(ii,:) = cell2mat(module_id_rgb(ind,3:5));
        mnames(ii) = module_id_rgb(ind,1);
    else
        mcolors(ii,:) = [0.7,0.7,0.7];
        mnames(ii) = {'ns'};
    end
end
br.CData = mcolors;        
 
% Axis properties
ax = gca;
ax.YLim = [-.5 13];
ax.YLabel.String = 'Mean Q*_c (z-score)';
ax.XTick = 1:size(data,2);
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 270;
ax.Title.String = 'Symbolic Community Significance';
ax.FontSize = 25;
grid on
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

% % Plot data points
% hold on
% for ii = 1:length(data_median)
%     for jj = 1:length(Qc_sig_zval_sort)
%       pl = plot(ii,Qc_sig_zval_sort(ii,jj),'.');
%       pl.MarkerEdgeColor = 'k';
%     end
% end
% hold off

%% ~~~~~~ NONSYMBOLIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Load data
load('workspace_Qc_MODULESfrom_BN202_Compare_Nonsymbolic_DATAfrom_BN202_Compare_Nonsymbolic.mat',...
     'Qc_sig_zval','Qc_sig_pval','Qc_sig_hval','Qc_perm','Qc_true','cond_data','cond_modules','modules_true');
 
%% Get significant modules based on p-value (since some of the permuted distributions are not normal,i.e. have long positive tail)
sig_modules = modules_true(median(Qc_sig_pval,2) < 0.01);

%% Sort the modules by mean z-score
[~,mz_sort] = sort(mean(Qc_sig_zval,2),'descend');
modules_true_sort = modules_true(1,mz_sort);
Qc_sig_zval_sort = Qc_sig_zval(mz_sort,:);

%% Plot selected distributions
% for ii = 1:size(Qc_perm,3)
%     d = squeeze(Qc_perm(13,:,ii));
%     t = Qc_true(13,ii);
%     figure;
%     plot_pairperm_distribution(t,d,'',1,'');
% end

%% Plot results
clear mcolors mnames std_err_mean
% Custom bar chart
data = Qc_sig_zval_sort';
data_mean = mean(data);
for ii = 1:size(data,2)
    std_err_mean(1,ii) = nanstd(data(:,ii))/sqrt(numel(~isnan(data(:,ii))));
end
figure('Position',[100,100,1400,800]);
br = bar(data_mean);
br.BarWidth = 0.9;
br.FaceColor = 'flat';
% Set the module colors and names from the defined cell array
for ii = 1:length(data_mean)
    if ismember(modules_true_sort(ii),cell2mat(module_id_rgb(:,2))) && ismember(modules_true_sort(ii),sig_modules)
        ind = find(cell2mat(module_id_rgb(:,2)) == modules_true_sort(ii));
        mcolors(ii,:) = cell2mat(module_id_rgb(ind,3:5));
        mnames(ii) = module_id_rgb(ind,1);
    else
        mcolors(ii,:) = [0.7,0.7,0.7];
        mnames(ii) = {'ns'};
    end
end
br.CData = mcolors;        
 
% Axis properties
ax = gca;
ax.YLim = [-.5 13];
ax.YLabel.String = 'Mean Q*_c (z-score)';
ax.XTick = 1:size(data,2);
ax.XTickLabel = mnames;
ax.XTickLabelRotation = 270;
ax.Title.String = 'Nonsymbolic Community Significance';
ax.FontSize = 25;
grid on
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

% % Plot data points
% hold on
% for ii = 1:length(data_median)
%     for jj = 1:length(Qc_sig_zval_sort)
%       pl = plot(ii,Qc_sig_zval_sort(ii,jj),'.');
%       pl.MarkerEdgeColor = 'k';
%     end
% end
% hold off