purge;

cd 'Z:\CR7T_Connectivity\MRI_proc'
fnames = subdir('CR_*_stimvars+rt+acc.mat');

for ii = 1:numel(fnames)
    load(fnames(ii).name,'compare_A','compare_B');
    format = [compare_A.format;         compare_B.format];
    acc =    [compare_A.Stimuli0x2EACC; compare_B.Stimuli0x2EACC];
    rt =     [compare_A.Stimuli0x2ERT;  compare_B.Stimuli0x2ERT];
    resp =   [compare_A.Stimuli0x2ERESP;compare_B.Stimuli0x2ERESP];
    
    % Restrict to valid trials, i.e. where there was no tick mark from
    % scanner
    valid = resp ~= '{`}'; valid = valid(:,1);
    valid_S = format(valid) == 1; % Symbolic
    valid_N = format(valid) == 2; % Nonsymbolic
    valid_S_corr = format(valid) == 1 & acc(valid) == 1; % Symbolic + Correct
    valid_N_corr = format(valid) == 2 & acc(valid) == 1; % Nonsymbolic + Correct
    
    acc = acc(valid);
    rt = rt(valid);
    
    % Mean accuracy
    allsubs_acc_mean_S(ii,1) = mean(acc(valid_S));
    allsubs_acc_mean_N(ii,1) = mean(acc(valid_N));
    
    % Mean RT in all trials
    allsubs_rt_mean_S(ii,1) = mean(rt(valid_S));
    allsubs_rt_mean_N(ii,1) = mean(rt(valid_N));
    
    % Mean RT in correct trials
    allsubs_rt_corr_mean_S(ii,1) = mean(rt(valid_S_corr));
    allsubs_rt_corr_mean_N(ii,1) = mean(rt(valid_N_corr));
    
    % Total trials valid
    allsubs_trials_valid_S(ii,1) = sum(valid_S);
    allsubs_trials_valid_N(ii,1) = sum(valid_N);
    
end

%% Paired Ttests
remove_subs = [5,6,8,11,12,28,32];
reduced_list = logical(ones(length(fnames),1));
reduced_list(remove_subs) = 0;

[h,p_rt,ci,stat_rt] =   ttest(allsubs_rt_corr_mean_S(reduced_list), allsubs_rt_corr_mean_N(reduced_list))
[h,p_acc,ci,stat_acc] = ttest(allsubs_acc_mean_S(reduced_list), allsubs_acc_mean_N(reduced_list))


[h,p_rt,ci,stat_rt] =   bf.ttest(allsubs_rt_corr_mean_S(reduced_list), allsubs_rt_corr_mean_N(reduced_list))
[h,p_acc,ci,stat_acc] = bf.ttest(allsubs_acc_mean_S(reduced_list), allsubs_acc_mean_N(reduced_list))


%% Construct values for Table 1
data = [allsubs_acc_mean_S(reduced_list),...
       allsubs_acc_mean_N(reduced_list),...
       allsubs_rt_corr_mean_S(reduced_list),...
       allsubs_rt_corr_mean_N(reduced_list)]; 
for ii = 1:size(data,2)
    d = data(:,ii);
    table1(ii,1) = mean(d);
    table1(ii,2) = std(d);
    table1(ii,3) = min(d);
    table1(ii,4) = max(d);
end

%% Plotting
% Connected boxplot - Accuracy
data1 = allsubs_acc_mean_S(reduced_list);
data2 = allsubs_acc_mean_N(reduced_list);
figure
boxplot([data1,data2]);
hold on
for ii=1:length(data1)
  plot([1,2],[data1(ii),data2(ii)],'-o',...
       'Color', [0.5,0.5,0.5],...
       'MarkerFaceColor',[1,0.5,0.5],...
       'MarkerEdgeColor',[0.5,0.5,0.5]);
end
hold off
ax = gca;
ax.XTickLabel = {'Symbolic';'Nonsymbolic'};
title('Mean Accuracy - Symbolic vs Nonsymbolic Trials');
ylabel('Accuracy');
set(gca,'FontSize',10);
grid on

% Connected boxplot - Accuracy
data1 = allsubs_rt_corr_mean_S(reduced_list);
data2 = allsubs_rt_corr_mean_N(reduced_list);
figure
boxplot([data1,data2]);
hold on
for ii=1:length(data1)
  plot([1,2],[data1(ii),data2(ii)],'-o',...
       'Color', [0.5,0.5,0.5],...
       'MarkerFaceColor',[1,0.5,0.5],...
       'MarkerEdgeColor',[0.5,0.5,0.5]);
end
hold off
ax = gca;
ax.XTickLabel = {'Symbolic';'Nonsymbolic'};
title('Mean Reaction Time - Symbolic vs Nonsymbolic Trials');
ylabel('Reaction Time (ms)');
set(gca,'FontSize',10);
grid on
    