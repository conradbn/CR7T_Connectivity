%% Modularity Q* differences between formats
purge

load('workspace_consensus_clustering_gamma_2pt45.mat','Q_all','Pa_all','Pa_group','zscore_mats')

Q_sym = squeeze(Q_all.BN202_Compare_Symbolic);
Q_nonsym = squeeze(Q_all.BN202_Compare_Nonsymbolic);
Pa_sym = squeeze(Pa_all.BN202_Compare_Symbolic);
Pa_nonsym = squeeze(Pa_all.BN202_Compare_Nonsymbolic);

% Get max Q* values
for ii = 1:size(Q_sym,2)
    %[h(ii),p(ii),ci(ii),stats(ii)] = ttest2(Q_sym(:,ii),Q_nonsym(:,ii));
    [h(ii),p(ii)] = ttest2(Q_sym(:,ii),Q_nonsym(:,ii));
    Q_nonsym_med(ii,1) = median(Q_nonsym(:,ii));
    Q_sym_med(ii,1) = median(Q_sym(:,ii));
    Q_nonsym_max(ii,1) = max(Q_nonsym(:,ii));
    Q_sym_max(ii,1) = max(Q_sym(:,ii));
end

% Connected boxplot - max Qstar
figure
boxplot([Q_sym_max,Q_nonsym_max]);
hold on
for ii=1:length(Q_sym_max)
  plot([1,2],[Q_sym_max(ii),Q_nonsym_max(ii)],'-o',...
       'Color', [0.5,0.5,0.5],...
       'MarkerFaceColor',[1,0.5,0.5],...
       'MarkerEdgeColor',[0.5,0.5,0.5]);
end
hold off
ax = gca;
ax.XTickLabel = {'Symbolic';'Nonsymbolic'};
ax.Title.String = 'Global Modularity in Subject-Level Connectivity Matrices';
ax.YLabel.String = 'Max Q* (modularity index)';
ax.FontSize = 14; 
grid on

% Paired Ttest
[h,p,ci,stats] = ttest(Q_sym_max,Q_nonsym_max)
bfac = bf.ttest(Q_sym_max,Q_nonsym_max)
[p,h,stat] = signrank(Q_sym_max,Q_nonsym_max)
















    