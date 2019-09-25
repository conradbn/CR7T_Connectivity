%% Run on individual subject matrices
purge
cd Z:\CR7T_Connectivity\matfiles
load('workspace_zscore_mats_33subs.mat','zscore_mats');
mats_struct = zscore_mats;
%mats = mats_struct.BN199_Compare_All;
gammas = 2.45;
tau = 0.5;
niters = 1000;

% Run the function
[C_all,M_all,Q_all,Pa_all,Pa_group,C_group] = ...
    CR7T_louvain_with_consensus_clustering_group_partition(mats_struct,gammas,tau,niters);

% Save workspace
save('workspace_consensus_clustering_gamma_2pt45.mat');