%% Calculate flexibility over gamma sweep
% Flexibility calculated using separate group level consensus partitions of 
% Symbolic and Nonsymbolic Zscore matrices

purge
load('workspace_modularity_gamma_sweep_33subs_Symbolic.mat','C_group_iter','gammas','iters_of_each_gamma');
Cs = C_group_iter;
load('workspace_modularity_gamma_sweep_33subs_Nonsymbolic.mat','C_group_iter','gammas','iters_of_each_gamma');
Cn = C_group_iter;

%% Get set of all unique pairs of partitions (note nchoosek does not give the repetitions)
pairs = nchoosek(1:size(Cs,3),2);
% Add the repetitions of each position (i.e. 1 and 1)
reps = repmat(1:size(Cs,3),2,1)';
pairs = [pairs;reps];

%% Loop through each level of gamma
for ii = 1:numel(gammas)
    disp(num2str(gammas(ii)));
    parfor jj = 1:length(pairs)
        % Set temporary values
        s = Cs(:,ii,pairs(jj));
        n = Cn(:,ii,pairs(jj));

        % Relabel so as to match community assignments as close as possible
        s2 = pair_labeling(n,s);
        % Combine partition vectors into format for flexibility function
        Cpair = [s2,n]';
        flex1 = flexibility(Cpair,'cat');
        flex(:,ii,jj) = flex1;
        
        % Alternative flexibility calculation
        % If singletons in either condition give value of 1
        flex2 = flex1;
        singleton = find(histc(n,unique(n))==1);
        flex2(ismember(n,singleton)) = 1;
        singleton = find(histc(s2,unique(s2))==1);
        flex2(ismember(s2,singleton)) = 1;
        flex_alt(:,ii,jj) = flex2;
    end
end

%% Count number of unique modules at each gamma
for ii = 1:numel(gammas)
    disp(num2str(gammas(ii)));
    for jj = 1:iters_of_each_gamma
        % Count unique values
        s = Cs(:,ii,jj);
        n = Cn(:,ii,jj);
        Cs_count(ii,jj) = numel(unique(s));
        Cn_count(ii,jj) = numel(unique(n));
        % Exclude singleton communities (regions by themself)
        singleton = find(histc(s,unique(s))==1);
        Cs_count_rm_sngltn(ii,jj) = numel(unique(s)) - numel(singleton);
        Cs_count_sngltn(ii,jj) = numel(singleton); % Capture the # of singletons
        
        singleton = find(histc(n,unique(n))==1);
        Cn_count_rm_sngltn(ii,jj) = numel(unique(n)) - numel(singleton);
        Cn_count_sngltn(ii,jj) = numel(singleton); % Capture the # of singletons

    end
end


%% Get the mena and std deviation (variability) of flexibility across pairwise comparisons at each gamma
flex_mean = squeeze(mean(flex,1));
flex_std = squeeze(std(flex,1));
% Get the mean variability at each gamma
flex_std_mean = mean(flex_std,2);
flex_mean_mean = mean(flex_mean,2);
% Get the gamma (index) at which the variability mean is maximal
flex_std_mean_max_index = find(flex_std_mean == max(flex_std_mean));

flexB_std = squeeze(std(flex_alt,1));
% Get the mean variability at each gamma
flexB_std_mean = mean(flexB_std,2);
% Get the gamma (index) at which the variability mean is maximal
flexB_std_mean_max_index = find(flexB_std_mean == max(flexB_std_mean));

%% Save workspace
%save('workspace_flexibility_sym_nonsym_gamma_sweep.mat');

%% Plotting - Variability of Flexibility over gamma sweep
% Default boxplot of standard deviations
figure
notBoxPlot(flex_std');
title('Variability of Flexibility Across Nodes - Gamma Sweep');
grid on
% Plot of mean variability over gamma
figure
plot(flex_std_mean);
title('Mean Variability of Flexibility Across Nodes - Gamma Sweep');
grid on

figure
plot(flexB_std_mean);
title('Mean Variability of Flexibility (Alt) Across Nodes - Gamma Sweep');
grid on

%% Custom boxplots - Variability of Flexibility
import iosr.statistics.*
figure('Position',[100,100,1200,550])
% Plot red line at peak
pl = line([flex_std_mean_max_index flex_std_mean_max_index], [0 0.5]);
pl.LineWidth = 2;
pl.LineStyle = '--';
pl.Color = 'r';
hold on
% Boxplot customization
bp = boxPlot(flex_std');
bp.boxColor = [.7 .7 .7];
bp.lineWidth = 2;
bp.lineColor = [.7 .7 .7];
bp.medianColor = [.7 .7 .7];
bp.showMean = true;
bp.meanMarker = '.';
bp.meanColor = 'k';
bp.meanSize = 20;
bp.showOutliers = false;
% bp.symbolMarker = '.'
% bp.symbolColor = [.8 .8 .8]
grid on
ax = gca;
ax.FontSize = 17;
ax.YLim = [0.1 0.5];
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
ax.YLabel.String = 'Standard Deviation of Flexibility';


%% Custom boxplots - Number of modules
figure('Position',[100,100,550,550]);
pl = plot(mean(Cs_count'));
pl.LineWidth = 3;
pl.Color = 'r';
hold on; grid on;
pl = plot(mean(Cn_count'));
pl.LineWidth = 3;
pl.Color = 'b';
ax = gca;
ax.FontSize = 17;
ax.XTick = [0 10 20 30 40 50 60 70 80 90 100];
ax.XTickLabel = [{'0'} num2cell(gammas(ax.XTick(2:end)))];
ax.XTickLabelRotation = 270;
ax.YLabel.String = 'Number of modules';
ax.XLabel.String = 'Gamma Value';
ax.Title.String = 'Number of Modules';
legend('Symbolic','Nonsymbolic','Location','southeast');

figure('Position',[100,100,550,550]);
pl = plot(mean(Cs_count_rm_sngltn'));
pl.LineWidth = 3;
pl.Color = 'r';
hold on; grid on;
pl = plot(mean(Cn_count_rm_sngltn'));
pl.LineWidth = 3;
pl.Color = 'b';
ax = gca;
ax.FontSize = 17;
ax.XTick = [0 10 20 30 40 50 60 70 80 90 100];
ax.XTickLabel = [{'0'} num2cell(gammas(ax.XTick(2:end)))];
ax.XTickLabelRotation = 270;
ax.YLabel.String = 'Number of modules';
ax.XLabel.String = 'Gamma Value';
ax.Title.String = 'Number of Non-Singleton Modules';
legend('Symbolic','Nonsymbolic','Location','southeast');

figure('Position',[100,100,550,550]);
pl = plot(mean(Cs_count_sngltn'));
pl.LineWidth = 3;
pl.Color = 'r';
hold on; grid on;
pl = plot(mean(Cn_count_sngltn'));
pl.LineWidth = 3;
pl.Color = 'b';
ax = gca;
ax.FontSize = 17;
ax.XTick = [0 10 20 30 40 50 60 70 80 90 100];
ax.XTickLabel = [{'0'} num2cell(gammas(ax.XTick(2:end)))];
ax.XTickLabelRotation = 270;
ax.YLabel.String = 'Number of modules';
ax.XLabel.String = 'Gamma Value';
ax.Title.String = 'Number of Singletons';
legend('Symbolic','Nonsymbolic','Location','southeast');


%% Custom boxplots - Variability of Flexibility (alternative calculation)
% import iosr.statistics.*
% figure('Position',[100,100,1400,500])
% % Boxplot customization
% bp = boxPlot(flexB_std');
% bp.boxColor = [.7 .7 .7];
% bp.lineWidth = 2;
% bp.lineColor = [.7 .7 .7];
% bp.medianColor = [.7 .7 .7];
% bp.showMean = true;
% bp.meanMarker = '.';
% bp.meanColor = 'k';
% bp.meanSize = 20;
% bp.showOutliers = false;
% % bp.symbolMarker = '.'
% % bp.symbolColor = [.8 .8 .8]
% grid on
% ax = gca;
% ax.YLim = [0 0.55]
% ax.FontSize = 12;
% %ax.GridAlpha = .3;
% % Title
% ax.Title.String = 'Variability of Flexibility (Alt) Across Nodes - Gamma Sweep';
% ax.Title.Interpreter = 'none';
% % X Ticks
% ax.XTick = 1:numel(gammas);
% gcells = num2cell(gammas);
% ax.XTickLabel = gcells(1:end);
% ax.XTickLabelRotation = 270;
% ax.XAxis.FontSize = 7;
% % Axis labels
% xL = xlabel('Gamma Value');
% xL.FontSize = 12;
% ax.YLabel.String = 'Standard Deviation of Flexibility';
% 
