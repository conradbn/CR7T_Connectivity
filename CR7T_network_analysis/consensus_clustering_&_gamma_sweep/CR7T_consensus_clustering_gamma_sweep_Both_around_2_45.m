%% Gamma sweep - Consensus clustering
% Run full group consensus over mutliple iterations at each gamma

% NOTE - This script is for the gamma sweep around at +/- 1 around gamma =
% 2.45. For the VAN/DAN/FPN robustness check in Supplementary Figure 2.
% 1000 iterations used at each gamma step to mimic that used for the
% primary analysis at gamma = 2.45.

%% NONSYMBOLIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
purge
% Load data
load('workspace_zscore_mats.mat','zscore_mats');
% Get mats for the gamma sweep analysis
mats_struct.BN202_Compare_Nonsymbolic = zscore_mats.BN202_Compare_Nonsymbolic;

% Exclude subjects
subs_exclude = [5,6,8,11,12,28,32];
mats_struct.BN202_Compare_Nonsymbolic(:,:,subs_exclude) = [];

% Set parameters
gammas = 1.45:0.05:3.45; % Gamma range
tau = 0.5; % Tau parameter for consensus clustering funciton
niters = 1000;%10; % Iterations of Louvain algorithm for each subject matrix
iters_of_each_gamma = 1; % Iterations of full group consensus at each gamma value

%% Run all iterations, using parallelization
%parfor_progress(iters_of_each_gamma);
for ii = 1:iters_of_each_gamma
    % Run the function (return only the final group consensus partitions, C_group)
    [~,~,~,~,~,C_group] = ...
        CR7T_louvain_with_consensus_clustering_group_partition(mats_struct,gammas,tau,niters);
     C_group_iter(:,:,ii) = C_group.BN202_Compare_Nonsymbolic;
     parfor_progress;
end
%parfor_progress(0);

% Save values (large file, need to use v7.3
%save('Z:\CR7T_Connectivity\matfiles\workspace_modularity_gamma_sweep_33subs_Nonsymbolic.mat','-v7.3');
%save('Z:\CR7T_Connectivity\matfiles\workspace_modularity_gamma_sweep_33subs_Nonsymbolic_100xLouvain.mat','-v7.3');
save('Z:\CR7T_Connectivity\matfiles\workspace_modularity_gamma_sweep_33subs_Nonsymbolic_1000xLouvain_1xGamma_2pt45.mat','-v7.3');

%% SYMBOLIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
purge
% Load data
load('workspace_zscore_mats.mat','zscore_mats');
% Get mats for the gamma sweep analysis
mats_struct.BN202_Compare_Symbolic = zscore_mats.BN202_Compare_Symbolic;

% Exclude subjects
subs_exclude = [5,6,8,11,12,28,32];
mats_struct.BN202_Compare_Symbolic(:,:,subs_exclude) = [];

% Set parameters
gammas = 1.45:0.05:3.45; % Gamma range
tau = 0.5; % Tau parameter for consensus clustering funciton
niters = 1000;%10; % Iterations of Louvain algorithm for each subject matrix
iters_of_each_gamma = 1; % Iterations of full group consensus at each gamma value

%% Run all iterations, using parallelization
%parfor_progress(iters_of_each_gamma);
for ii = 1:iters_of_each_gamma
    % Run the function (return only the final group consensus partitions, C_group)
    [~,~,~,~,~,C_group] = ...
        CR7T_louvain_with_consensus_clustering_group_partition(mats_struct,gammas,tau,niters);
     C_group_iter(:,:,ii) = C_group.BN202_Compare_Symbolic;
     parfor_progress;
end
%parfor_progress(0);

% Save values (large file, need to use v7.3
%save('Z:\CR7T_Connectivity\matfiles\workspace_modularity_gamma_sweep_33subs_Symbolic.mat','-v7.3');
%save('Z:\CR7T_Connectivity\matfiles\workspace_modularity_gamma_sweep_33subs_Symbolic_100xLouvain.mat','-v7.3');
save('Z:\CR7T_Connectivity\matfiles\workspace_modularity_gamma_sweep_33subs_Symbolic_1000xLouvain_1xGamma_2pt45.mat','-v7.3');
